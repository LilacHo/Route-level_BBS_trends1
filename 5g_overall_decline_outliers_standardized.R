# 5g_overall_decline_outliers_standardized.R
#
# Fixes the "common species problem" found in 5f_overall_decline_outliers.R:
# once a route's n_species pool is most of the local bird community (because
# the species filter now passes 62% of species, and the route filter accepts
# both Decreasing-suitability and Stable categories), a route's mean_trend is
# really "is the average local bird doing well here" -- and the "top species"
# behind a cluster just lists whichever generalists (American Robin, Barn
# Swallow, ...) happen to live in that region, present at nearly every route,
# contributing an unremarkable trend value each time. Common presence got
# confused with driving the signal.
#
# The fix (approach #2 of the three discussed): before averaging a route's
# species trends together, standardize each species-route trend against that
# SAME SPECIES' OWN across-route baseline (its mean and SD across all of its
# converged routes, any SDM category) --
#
#   z(species, route) = (trend(species, route) - mean_trend(species)) / sd_trend(species)
#
# A generalist with an unremarkable trend everywhere now contributes close to
# zero to a route's score no matter how often it shows up, because this route
# isn't unusual FOR THAT SPECIES. A species doing surprisingly well or badly
# specifically at this route (relative to how it normally does) contributes a
# lot. The route-level aggregate (mean_z) now means "are species here doing
# unusually well/badly relative to their own normal range," not "is the
# average bird doing well here" -- and the species tables rank by |z| (how
# anomalous that species is AT that route), not by occurrence count, so a
# species only shows up as "driving" a route if it's actually behaving
# unusually there.
#
# Steps 1 (species filter) and 2 (route filter) are otherwise identical to
# 5f: a species counts as an overall decliner only using its credible routes
# (own 90% CI excludes zero), and qualifying routes are Decreasing-suitability
# OR Stable (Increasing-suitability excluded). Invasive species excluded as
# in 5b/5e/5f. Outlier classification is still the plain percentile cutoff
# from 5f (no Gi*, no spatial smoothing, no multiple-testing correction).
#
# Reads:  output/species_routes_covariates/per_species_sdm/*_route_trends_sdm.csv
# Writes: output/species_routes_covariates/hotspots/overall_decline_outliers_zscore/<model_tag>_species_baseline.csv
#         output/species_routes_covariates/hotspots/overall_decline_outliers_zscore/<run_label>_<model_tag>_<scenario>_noinvasive_route_outliers.csv
#         output/species_routes_covariates/hotspots/overall_decline_outliers_zscore/<run_label>_<model_tag>_<scenario>_noinvasive_{increasing,decreasing}_species.csv
#         output/species_routes_covariates/hotspots/overall_decline_outliers_zscore/plots/<run_label>_<model_tag>_<scenario>_noinvasive_outlier_map.png

library(here)
library(tidyverse)
library(sf)

here::i_am("5g_overall_decline_outliers_standardized.R")

# Settings -------------------------------------------------------------------
bird_group <- NA   # NA pools every group; set to one Group value to restrict.

model_tags <- c("base", "anthro")
scenarios  <- c("rcp45", "rcp85")

require_route_converged <- TRUE

# Same invasive-species exclusion as 5b/5e/5f. Set to character(0) to disable.
exclude_species_codes <- {
  traits <- read.csv(here::here("data",
    "spp_names_codes_group_aou_taxa_redlist_diet_migration_invasive.csv"))
  traits$Code[traits$invasive == "Widespread invasive"]
}

# Step 1 (species filter, same as 5f): a species' overall trend is the mean
# of ONLY its credible routes (own 90% CI excludes zero), across all
# converged routes of any SDM category. Must be < 0.
min_credible_routes_overall <- 1

# Step 2 (route filter, same as 5f): SDM categories a species-route pair
# must fall into to be considered.
include_categories <- c("Decreasing suitability", "Stable")

# NEW -- species baseline for standardization: a species needs at least this
# many converged routes (any SDM category, same set species_baseline is
# built from) before its own mean/SD is trusted enough to standardize
# against. Too few routes gives an unstable SD (or SD = 0), which would
# blow up that species' z-scores -- those species are dropped from the
# z-score step entirely (but still count toward Step 1's species filter).
min_routes_for_baseline <- 5

# Percentile cutoffs for the outlier classification (two-tailed, nested),
# applied to the route-level mean_z instead of 5f's mean_trend.
outlier_probs <- c(0.01, 0.05, 0.10)

in_dir   <- here::here("output", "species_routes_covariates", "per_species_sdm")
out_dir  <- here::here("output", "species_routes_covariates", "hotspots", "overall_decline_outliers_zscore")
plot_dir <- file.path(out_dir, "plots")
if (!dir.exists(out_dir))  dir.create(out_dir,  recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Read every SDM CSV (same input as 5b/5d/5e/5f) ------------------------------
sdm_files <- list.files(in_dir, pattern = "_route_trends_sdm\\.csv$", full.names = TRUE)
if (length(sdm_files) == 0) {
  stop("No SDM CSVs found in ", in_dir, " — run 3c_add_SDM_covariates.R first.")
}

cat("\nCombining all species SDM files...\n")
all_sdm_raw_unfiltered <- sdm_files %>% map_dfr(read.csv)

if (is.na(bird_group)) {
  all_sdm_raw <- all_sdm_raw_unfiltered
  run_label   <- "all_groups"
} else {
  all_sdm_raw <- all_sdm_raw_unfiltered %>% filter(group == bird_group)
  run_label   <- bird_group
}

if (require_route_converged) {
  n_before <- nrow(all_sdm_raw)
  all_sdm_raw <- all_sdm_raw %>% filter(route_converged == TRUE)
  cat("require_route_converged = TRUE -> dropped", n_before - nrow(all_sdm_raw),
      "of", n_before, "row(s).\n")
}

if (length(exclude_species_codes) > 0) {
  n_before <- nrow(all_sdm_raw)
  excluded_present <- intersect(exclude_species_codes, unique(all_sdm_raw$species_code))
  all_sdm_raw <- all_sdm_raw %>% filter(!species_code %in% exclude_species_codes)
  cat("exclude_species_codes -> dropped", n_before - nrow(all_sdm_raw), "of", n_before,
      "row(s) (", paste(excluded_present, collapse = ", "), ")\n")
}

if (nrow(all_sdm_raw) == 0) {
  stop("No rows left after filtering (bird_group = '", bird_group,
       "', require_route_converged = ", require_route_converged, ").")
}

# Decreasing suitability = 1,2,3 | Stable = 4 | Increasing suitability = 5,6,7
group_category <- function(x) {
  dplyr::case_when(
    x %in% c(1, 2, 3) ~ "Decreasing suitability",
    x == 4            ~ "Stable",
    x %in% c(5, 6, 7) ~ "Increasing suitability",
    TRUE              ~ NA_character_
  )
}

# US-states background map (same as 5_hotspots.R / 5b / 5e / 5f).
data("us_states", package = "spData")
us_states <- st_transform(us_states, 4326)
us_bbox <- st_bbox(us_states)
map_xlim <- c(us_bbox["xmin"] - 1, us_bbox["xmax"] + 1)
map_ylim <- c(us_bbox["ymin"] - 1, us_bbox["ymax"] + 1)

# Same plain-percentile classifier as 5f, just applied to mean_z.
classify_outliers <- function(x, probs = outlier_probs) {
  q_lo <- quantile(x, probs = probs, na.rm = TRUE, type = 7)
  q_hi <- quantile(x, probs = 1 - probs, na.rm = TRUE, type = 7)
  out <- dplyr::case_when(
    x <= q_lo[["1%"]]  ~ "Most Decreasing (99%)",
    x <= q_lo[["5%"]]  ~ "Most Decreasing (95%)",
    x <= q_lo[["10%"]] ~ "Most Decreasing (90%)",
    x >= q_hi[["99%"]] ~ "Most Increasing (99%)",
    x >= q_hi[["95%"]] ~ "Most Increasing (95%)",
    x >= q_hi[["90%"]] ~ "Most Increasing (90%)",
    TRUE               ~ "Not Outlier"
  )
  factor(out, levels = c("Most Decreasing (99%)", "Most Decreasing (95%)", "Most Decreasing (90%)",
                         "Not Outlier",
                         "Most Increasing (90%)", "Most Increasing (95%)", "Most Increasing (99%)"))
}

outlier_fill_colors <- c(
  "Most Decreasing (99%)" = "#2166AC", "Most Decreasing (95%)" = "#67A9CF", "Most Decreasing (90%)" = "#D1E5F0",
  "Not Outlier"            = "grey85",
  "Most Increasing (90%)"  = "#FDDBC7", "Most Increasing (95%)" = "#EF8A62", "Most Increasing (99%)" = "#B2182B"
)

make_outlier_map <- function(df, title, subtitle, file_path) {
  p <- ggplot() +
    geom_sf(data = us_states, fill = "grey97", color = "grey70", linewidth = 0.3) +
    geom_point(data = df, aes(x = longitude, y = latitude, color = outlier_class),
               size = 1.6, alpha = 0.85) +
    geom_point(data = df %>% filter(outlier_class != "Not Outlier"),
               aes(x = longitude, y = latitude),
               shape = 1, size = 2.8, stroke = 0.5, color = "black") +
    scale_color_manual(values = outlier_fill_colors, name = "Outlier\nclass", drop = FALSE) +
    coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
    labs(title = title, subtitle = subtitle,
         caption = "Value is each route's mean standardized anomaly (z-score vs. each species' own baseline), not raw trend -- no spatial smoothing.",
         x = NULL, y = NULL) +
    theme_minimal() +
    theme(plot.title    = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 10),
          plot.caption  = element_text(size = 8, color = "grey40"),
          axis.text     = element_text(size = 8))
  ggsave(file_path, p, width = 9, height = 6.5, dpi = 150)
  cat("Saved plot:", basename(file_path), "\n")
}

# Main loop: model_tag x scenario ---------------------------------------------
for (model_tag in model_tags) {
  target_sdm <- all_sdm_raw %>% filter(model == model_tag)
  if (nrow(target_sdm) == 0) {
    message("Skipping model_tag = '", model_tag, "' -- no rows found.")
    next
  }
  cat("\n=== model_tag:", model_tag, "===\n")

  # Step 1 -- species filter: overall trend from CREDIBLE routes only (any
  # SDM category), must be < 0. Identical to 5f.
  credible_rows <- target_sdm %>% filter(!is.na(trend), trend_lci > 0 | trend_uci < 0)
  n_before <- target_sdm %>% filter(!is.na(trend)) %>% nrow()
  cat("Credible routes (own 90% CI excludes zero):", nrow(credible_rows), "of", n_before, "\n")

  species_overall <- credible_rows %>%
    group_by(species, species_code) %>%
    summarise(overall_trend = mean(trend), n_credible_routes = n(), .groups = "drop") %>%
    filter(n_credible_routes >= min_credible_routes_overall)

  declining_species <- species_overall %>% filter(overall_trend < 0)
  cat("Species with a credible overall decreasing mean trend:", nrow(declining_species),
      "of", nrow(species_overall), "\n")

  # NEW -- species baseline for standardization: mean/SD of trend across ALL
  # of that species' converged routes (any SDM category) -- the same scope
  # Step 1 uses, but from every row (not just credible ones), since we want
  # this species' normal RANGE of route-level trends, not just its confident
  # ones.
  species_baseline <- target_sdm %>%
    filter(!is.na(trend)) %>%
    group_by(species, species_code) %>%
    summarise(mean_sp = mean(trend), sd_sp = sd(trend), n_sp = n(), .groups = "drop") %>%
    filter(n_sp >= min_routes_for_baseline, !is.na(sd_sp), sd_sp > 0)

  cat("Species with a usable baseline (>=", min_routes_for_baseline, "routes, SD > 0):",
      nrow(species_baseline), "\n")

  write.csv(species_baseline,
            file.path(out_dir, paste0(run_label, "_", model_tag, "_species_baseline.csv")),
            row.names = FALSE)

  for (scenario in scenarios) {

    analysis <- target_sdm %>%
      filter(!is.na(.data[[scenario]]), !is.na(trend)) %>%
      transmute(species, species_code, group, route, latitude, longitude, trend,
                category = group_category(.data[[scenario]]))

    cat("\n----- model_tag:", model_tag, "| scenario:", scenario, "-----\n")

    # Step 2 -- route filter: Decreasing-suitability OR Stable routes for
    # overall-declining species only.
    qual_rows <- analysis %>%
      filter(category %in% include_categories, species_code %in% declining_species$species_code)

    if (nrow(qual_rows) == 0) {
      message("  No qualifying rows -- skipping.")
      next
    }
    cat("  Decreasing-suitability/Stable rows for overall-declining species:", nrow(qual_rows), "\n")

    # Standardize: join each row to its species' baseline and compute z.
    # Rows for species without a usable baseline (too few routes, or SD = 0)
    # are dropped here -- they can't be standardized, so they can't
    # contribute to mean_z (this is the price of fixing the common-species
    # problem: a species needs enough of its own route history to tell
    # whether it's doing something unusual here).
    n_before_z <- nrow(qual_rows)
    qual_rows <- qual_rows %>%
      inner_join(species_baseline %>% select(species_code, mean_sp, sd_sp, n_sp), by = "species_code") %>%
      mutate(z = (trend - mean_sp) / sd_sp)
    cat("  Rows with a usable species baseline:", nrow(qual_rows), "of", n_before_z, "\n")

    if (nrow(qual_rows) == 0) {
      message("  No rows with a usable baseline -- skipping.")
      next
    }

    route_summary <- qual_rows %>%
      group_by(route) %>%
      summarise(latitude = first(latitude), longitude = first(longitude),
                n_species = n(), mean_trend = mean(trend), mean_z = mean(z), .groups = "drop")

    cat("  Eligible routes:", nrow(route_summary), "\n")

    route_summary$outlier_class <- classify_outliers(route_summary$mean_z)

    n_inc <- sum(grepl("^Most Increasing", route_summary$outlier_class))
    n_dec <- sum(grepl("^Most Decreasing", route_summary$outlier_class))
    cat("  Most Increasing outliers (any tier):", n_inc, "\n")
    cat("  Most Decreasing outliers (any tier):", n_dec, "\n")

    file_prefix <- paste0(run_label, "_", model_tag, "_", scenario,
                          if (length(exclude_species_codes) > 0) "_noinvasive" else "")

    out_csv <- file.path(out_dir, paste0(file_prefix, "_route_outliers.csv"))
    write.csv(route_summary, out_csv, row.names = FALSE)
    cat("  Wrote:", basename(out_csv), "\n")

    # Per-route species detail behind the outlier routes, ranked by |z| (how
    # anomalous that species is AT that route) instead of by occurrence
    # count -- this is the actual fix for the common-species problem.
    increasing_species <- qual_rows %>%
      filter(route %in% route_summary$route[grepl("^Most Increasing", route_summary$outlier_class)]) %>%
      select(route, species, group, trend, z) %>%
      arrange(route, desc(z))
    write.csv(increasing_species, file.path(out_dir, paste0(file_prefix, "_increasing_species.csv")),
              row.names = FALSE)

    decreasing_species <- qual_rows %>%
      filter(route %in% route_summary$route[grepl("^Most Decreasing", route_summary$outlier_class)]) %>%
      select(route, species, group, trend, z) %>%
      arrange(route, z)
    write.csv(decreasing_species, file.path(out_dir, paste0(file_prefix, "_decreasing_species.csv")),
              row.names = FALSE)

    model_scenario_txt <- paste0(model_tag, ", ", toupper(sub("rcp", "RCP ", scenario)))
    make_outlier_map(
      route_summary,
      paste0("Standardized anomaly outliers, overall-declining species (", model_scenario_txt, ")"),
      paste0(nrow(route_summary), " eligible routes | ", n_inc, " most-increasing | ", n_dec, " most-decreasing outliers"),
      file.path(plot_dir, paste0(file_prefix, "_outlier_map.png"))
    )
  }
}

cat("\n=== Standardized outlier analysis complete (model_tags = ", paste(model_tags, collapse = ", "),
    ", scenarios = ", paste(scenarios, collapse = ", "), ", bird_group = ",
    if (is.na(bird_group)) "all_groups (NA)" else bird_group, ") ===\n", sep = "")
