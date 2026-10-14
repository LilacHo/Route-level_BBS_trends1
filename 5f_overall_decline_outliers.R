# 5f_overall_decline_outliers.R
#
# Restart of the 5e_overall_decline_hotspots.R idea, with two changes the
# user asked for after looking at 5e's results:
#
#   1. No more Getis-Ord Gi* / no spatial smoothing. Gi* asks "is this route
#      surrounded by other high/low routes" -- useful for finding coherent
#      REGIONS, but it also means a real standout route can get washed out
#      just because its neighbors are unremarkable, and conversely a
#      middling route can get flagged just for sitting next to one or two
#      extreme neighbors. This script instead ranks routes directly against
#      the WHOLE eligible set (no neighbors involved at all) and flags
#      outliers by empirical percentile -- the 90th/95th/99th most extreme
#      routes in each direction, same "90%/95%/99%" language 5b/5e used for
#      Gi* confidence classes, but now just a plain distribution cutoff.
#      There's no multiple-testing correction here (unlike Gi*'s FDR q-value)
#      because this isn't a hypothesis test -- it's a descriptive ranking.
#
#   2. The species filter is now always the stricter "credible decreaser":
#      a species' overall trend is the mean of ONLY its credible routes
#      (trend_lci > 0 or trend_uci < 0, any SDM category, across all its
#      converged routes), and that mean must be < 0. (5e tested this as an
#      optional toggle against a plain point-estimate version; this script
#      just uses the credible version directly, since 5e found it gave
#      nearly the same species pool anyway.)
#
#   3. The route filter now keeps BOTH SDM categories consistent with a
#      species getting worse or holding steady there -- "Decreasing
#      suitability" (1-3) AND "Stable" (4) -- rather than Decreasing
#      suitability alone. A route predicted Stable for a species that's
#      overall declining is still a meaningful place to ask "is this species
#      doing unusually well or badly here"; only "Increasing suitability"
#      (5-7) is excluded, since that's the one category where decline
#      wouldn't be surprising. The route's own trend value is NOT filtered
#      by its own credibility or sign -- it stays continuous (can be
#      increasing or decreasing), same as 5e.
#
# Same invasive-species exclusion as 5b/5e (this project's authoritative
# trait file, not a hand-picked guess).
#
# Reads:  output/species_routes_covariates/per_species_sdm/*_route_trends_sdm.csv
#         (written by 3c_add_SDM_covariates.R) -- same input as 5b/5d/5e.
# Writes: output/species_routes_covariates/hotspots/overall_decline_outliers/<model_tag>_species_overall_trend.csv
#         output/species_routes_covariates/hotspots/overall_decline_outliers/<run_label>_<model_tag>_<scenario>_noinvasive_route_outliers.csv
#         output/species_routes_covariates/hotspots/overall_decline_outliers/<run_label>_<model_tag>_<scenario>_noinvasive_{increasing,decreasing}_species.csv
#         output/species_routes_covariates/hotspots/overall_decline_outliers/plots/<run_label>_<model_tag>_<scenario>_noinvasive_outlier_map.png

library(here)
library(tidyverse)
library(sf)

here::i_am("5f_overall_decline_outliers.R")

# Settings -------------------------------------------------------------------
bird_group <- NA   # NA pools every group; set to one Group value to restrict.

model_tags <- c("base", "anthro")
scenarios  <- c("rcp45", "rcp85")

require_route_converged <- TRUE

# Same invasive-species exclusion as 5b/5e. Set to character(0) to disable.
exclude_species_codes <- {
  traits <- read.csv(here::here("data",
    "spp_names_codes_group_aou_taxa_redlist_diet_migration_invasive.csv"))
  traits$Code[traits$invasive == "Widespread invasive"]
}

# A species' overall trend is the mean of ONLY its credible routes (own 90%
# CI excludes zero), across all converged routes of any SDM category. That
# mean must be < 0 for the species to count as an overall decliner. See
# header comment point 2 for why this isn't a toggle here (5e tested both
# and they gave nearly the same species pool).
#
# Minimum number of credible routes a species needs before its overall trend
# is trusted at all.
min_credible_routes_overall <- 1

# Route filter: SDM categories to keep (a species' prediction at a route).
# "Decreasing suitability" and "Stable" both qualify -- see header comment
# point 3. "Increasing suitability" is always excluded.
include_categories <- c("Decreasing suitability", "Stable")

# Percentile cutoffs for the outlier classification (two-tailed, nested):
# 0.10/0.90 = 90% cutoff, 0.05/0.95 = 95%, 0.01/0.99 = 99%. Computed on the
# eligible routes' mean_trend directly -- no neighbors, no spatial weights.
outlier_probs <- c(0.01, 0.05, 0.10)

in_dir   <- here::here("output", "species_routes_covariates", "per_species_sdm")
out_dir  <- here::here("output", "species_routes_covariates", "hotspots", "overall_decline_outliers")
plot_dir <- file.path(out_dir, "plots")
if (!dir.exists(out_dir))  dir.create(out_dir,  recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Read every SDM CSV (same input as 5b/5d/5e) ---------------------------------
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
# -- kept in sync with 4c/5b/5d/5e.
group_category <- function(x) {
  dplyr::case_when(
    x %in% c(1, 2, 3) ~ "Decreasing suitability",
    x == 4            ~ "Stable",
    x %in% c(5, 6, 7) ~ "Increasing suitability",
    TRUE              ~ NA_character_
  )
}

# US-states background map (same as 5_hotspots.R / 5b / 5e).
data("us_states", package = "spData")
us_states <- st_transform(us_states, 4326)
us_bbox <- st_bbox(us_states)
map_xlim <- c(us_bbox["xmin"] - 1, us_bbox["xmax"] + 1)
map_ylim <- c(us_bbox["ymin"] - 1, us_bbox["ymax"] + 1)

# Classify each route's mean_trend against the eligible set's OWN empirical
# percentiles -- a plain distribution cutoff, no neighbors / no spatial
# weights, no multiple-testing correction (this is a ranking, not a test).
classify_outliers <- function(x, probs = outlier_probs) {
  q_lo <- quantile(x, probs = probs, na.rm = TRUE, type = 7)          # named 1%, 5%, 10%
  q_hi <- quantile(x, probs = 1 - probs, na.rm = TRUE, type = 7)      # named 99%, 95%, 90%
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
         caption = "Red = most increasing (doing unusually well); blue = most decreasing. Percentile cutoffs computed directly on this route set -- no spatial smoothing.",
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

  # Step 1 -- species filter: overall trend from CREDIBLE routes only (any
  # SDM category), must be < 0.
  credible_rows <- target_sdm %>% filter(!is.na(trend), trend_lci > 0 | trend_uci < 0)
  n_before <- target_sdm %>% filter(!is.na(trend)) %>% nrow()
  cat("\n=== model_tag:", model_tag, "===\n")
  cat("Credible routes (own 90% CI excludes zero):", nrow(credible_rows), "of", n_before, "\n")

  species_overall <- credible_rows %>%
    group_by(species, species_code) %>%
    summarise(overall_trend = mean(trend), n_credible_routes = n(), .groups = "drop") %>%
    filter(n_credible_routes >= min_credible_routes_overall)

  declining_species <- species_overall %>% filter(overall_trend < 0)
  cat("Species with a credible overall decreasing mean trend:", nrow(declining_species),
      "of", nrow(species_overall), "\n")

  write.csv(species_overall,
            file.path(out_dir, paste0(run_label, "_", model_tag, "_species_overall_trend.csv")),
            row.names = FALSE)

  for (scenario in scenarios) {

    analysis <- target_sdm %>%
      filter(!is.na(.data[[scenario]]), !is.na(trend)) %>%
      transmute(species, species_code, group, route, latitude, longitude, trend,
                category = group_category(.data[[scenario]]))

    cat("\n----- model_tag:", model_tag, "| scenario:", scenario, "-----\n")

    # Step 2 -- route filter: Decreasing-suitability OR Stable routes for
    # overall-declining species only. Route's own trend stays continuous
    # (not filtered by sign or credibility).
    qual_rows <- analysis %>%
      filter(category %in% include_categories, species_code %in% declining_species$species_code)

    if (nrow(qual_rows) == 0) {
      message("  No qualifying rows -- skipping.")
      next
    }
    cat("  Decreasing-suitability/Stable rows for overall-declining species:", nrow(qual_rows), "\n")

    route_summary <- qual_rows %>%
      group_by(route) %>%
      summarise(latitude = first(latitude), longitude = first(longitude),
                n_species = n(), mean_trend = mean(trend), .groups = "drop")

    cat("  Eligible routes:", nrow(route_summary), "\n")

    route_summary$outlier_class <- classify_outliers(route_summary$mean_trend)

    n_inc <- sum(grepl("^Most Increasing", route_summary$outlier_class))
    n_dec <- sum(grepl("^Most Decreasing", route_summary$outlier_class))
    cat("  Most Increasing outliers (any tier):", n_inc, "\n")
    cat("  Most Decreasing outliers (any tier):", n_dec, "\n")
    print(table(route_summary$outlier_class))

    file_prefix <- paste0(run_label, "_", model_tag, "_", scenario,
                          if (length(exclude_species_codes) > 0) "_noinvasive" else "")

    out_csv <- file.path(out_dir, paste0(file_prefix, "_route_outliers.csv"))
    write.csv(route_summary, out_csv, row.names = FALSE)
    cat("  Wrote:", basename(out_csv), "\n")

    # Per-route species detail behind the outlier routes in each direction.
    increasing_species <- qual_rows %>%
      filter(route %in% route_summary$route[grepl("^Most Increasing", route_summary$outlier_class)]) %>%
      count(route, species, group, sort = TRUE)
    write.csv(increasing_species, file.path(out_dir, paste0(file_prefix, "_increasing_species.csv")),
              row.names = FALSE)

    decreasing_species <- qual_rows %>%
      filter(route %in% route_summary$route[grepl("^Most Decreasing", route_summary$outlier_class)]) %>%
      count(route, species, group, sort = TRUE)
    write.csv(decreasing_species, file.path(out_dir, paste0(file_prefix, "_decreasing_species.csv")),
              row.names = FALSE)

    model_scenario_txt <- paste0(model_tag, ", ", toupper(sub("rcp", "RCP ", scenario)))
    make_outlier_map(
      route_summary,
      paste0("Overall-declining species, Decreasing-suitability/Stable routes (", model_scenario_txt, ")"),
      paste0(nrow(route_summary), " eligible routes | ", n_inc, " most-increasing | ", n_dec, " most-decreasing outliers"),
      file.path(plot_dir, paste0(file_prefix, "_outlier_map.png"))
    )
  }
}

cat("\n=== Overall-decline outlier analysis complete (model_tags = ", paste(model_tags, collapse = ", "),
    ", scenarios = ", paste(scenarios, collapse = ", "), ", bird_group = ",
    if (is.na(bird_group)) "all_groups (NA)" else bird_group, ") ===\n", sep = "")
