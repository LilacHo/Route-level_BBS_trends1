# 5e_overall_decline_hotspots.R
#
# A third framing of the refugia question, alongside 5b_refugia_hotspot_stats.R
# and 5d_species_voting_hotspots.R. Both of those start from a per-ROUTE
# criterion (is this species credibly increasing at this one route despite a
# Decreasing suitability prediction there?). This script starts from a per-SPECIES
# criterion instead:
#
#   1. Is the species an overall decliner? -- its mean trend across ALL of its
#      converged routes (every SDM category, not just Decreasing suitability -- the same
#      "overall" scope 6a_species_trend_distributions.R's dashed outline
#      uses), by point estimate only (no credible-interval requirement).
#   2. Restrict to that species' Decreasing-suitability routes only.
#   3. Unlike 5b/5d, do NOT then filter by the route's own trend sign -- a
#      route can be increasing or decreasing. Run Getis-Ord Gi* directly on
#      the route-level mean trend (continuous), among overall-declining
#      species' Decreasing-suitability routes.
#
# This Gi* run is two-sided by construction and both directions are
# meaningful here:
#   - Hot spots (z > 0): routes where overall-declining species are, on
#     average, doing LOCALLY BETTER than their neighbors despite being in
#     predicted-decreasing suitability territory -- the refugia candidates.
#   - Cold spots (z < 0): routes doing locally WORSE -- "compound decline"
#     zones where real population decline, predicted decreasing suitability, and local
#     underperformance all stack up.
#
# As with 5b, thousands of routes are tested at once, so Benjamini-Hochberg
# FDR-adjusted q-values (gi_q) are the basis for the headline counts, not the
# raw per-test 90/95/99% confidence classes (gi_class) alone.
#
# Reads:  output/species_routes_covariates/per_species_sdm/*_route_trends_sdm.csv
#         (written by 3c_add_SDM_covariates.R) -- same input as 5b/5d.
# Writes: output/species_routes_covariates/hotspots/overall_decline_stats/<model_tag>_<scenario>_species_overall_trend.csv
#         output/species_routes_covariates/hotspots/overall_decline_stats/<model_tag>_<scenario>_route_gi.csv
#         output/species_routes_covariates/hotspots/overall_decline_stats/plots/<model_tag>_<scenario>_gi_map.png

library(here)
library(tidyverse)
library(sf)
library(spdep)

here::i_am("5e_overall_decline_hotspots.R")

# Settings -------------------------------------------------------------------
bird_group <- NA   # NA pools every group; set to one Group value to restrict.

model_tags <- c("base", "anthro")
scenarios  <- c("rcp45", "rcp85")

require_route_converged <- TRUE

# Same invasive-species exclusion as 5b_refugia_hotspot_stats.R, applied
# before both the species-overall-trend step and the route-level step --
# see that script's header comment for why (Eurasian Collared-Dove alone
# behind most of the pooled magnitude hotspot before this was added). Pulled
# from this project's authoritative trait file, not a hand-picked guess. Set
# to character(0) to disable.
exclude_species_codes <- {
  traits <- read.csv(here::here("data",
    "spp_names_codes_group_aou_taxa_redlist_diet_migration_invasive.csv"))
  traits$Code[traits$invasive == "Widespread invasive"]
}

# FALSE reproduces the first test of this idea exactly: a species-route's own
# trend is used as-is (point estimate only), whatever its own credible
# interval looks like. TRUE restricts the Decreasing-suitability rows that
# feed the route-level mean (both the refugia and compound-decline
# directions) to rows whose OWN 90% credible interval excludes zero
# (trend_lci > 0 or trend_uci < 0) -- same "credible" convention 5_hotspots.R
# and 5b use, applied here to filter input rows rather than to define a
# route-level mismatch. Written to separately-named output (file_prefix gets
# a "_credible" suffix) so both versions sit side by side.
require_credible_route_trend <- FALSE

# FALSE reproduces the first test of this idea exactly: a species' overall
# trend is the plain mean of its route trends (point estimate only, any SDM
# category), whatever each route's own credible interval looks like. TRUE
# computes that same mean using only routes whose OWN 90% credible interval
# excludes zero (trend_lci > 0 or trend_uci < 0) -- a stricter definition of
# "overall decliner" that only counts a route's contribution to a species'
# trend when that route is itself a confident estimate. Independent of
# require_credible_route_trend above: that one filters which Decreasing suitability-
# predicted rows feed the route-level Gi* value; this one filters which rows
# feed the species-level classification step that happens before it. Output
# gets an "_overallcredible" suffix so this sits alongside the other variants.
require_credible_for_overall_trend <- FALSE

# Minimum number of converged routes (any SDM category) a species needs
# before its "overall trend" is trusted at all. 1 reproduces the first test
# of this idea exactly (no minimum); raise it to require more evidence behind
# the overall-decliner label.
min_routes_overall <- 1

# Number of nearest-neighbor routes for the Gi* spatial weights matrix (same
# rationale as 5b: routes aren't on a regular grid, so k-nearest-neighbors
# keeps neighborhood size comparable across variable route density).
k_neighbors <- 8

# Benjamini-Hochberg FDR threshold.
fdr_q <- 0.05

in_dir   <- here::here("output", "species_routes_covariates", "per_species_sdm")
out_dir  <- here::here("output", "species_routes_covariates", "hotspots", "overall_decline_stats")
plot_dir <- file.path(out_dir, "plots")
if (!dir.exists(out_dir))  dir.create(out_dir,  recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Read every SDM CSV (same input as 5b/5d) -----------------------------------
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

# Decreasing suitability = 1,2,3 | Stable = 4 | Increasing suitability = 5,6,7 -- kept in sync with
# 4c_statistical_analysis_and_visualization_covariates.R / 5b / 5d.
group_category <- function(x) {
  dplyr::case_when(
    x %in% c(1, 2, 3) ~ "Decreasing suitability",
    x == 4            ~ "Stable",
    x %in% c(5, 6, 7) ~ "Increasing suitability",
    TRUE              ~ NA_character_
  )
}

# US-states background map (same as 5_hotspots.R / 5b).
data("us_states", package = "spData")
us_states <- st_transform(us_states, 4326)
us_bbox <- st_bbox(us_states)
map_xlim <- c(us_bbox["xmin"] - 1, us_bbox["xmax"] + 1)
map_ylim <- c(us_bbox["ymin"] - 1, us_bbox["ymax"] + 1)

# Getis-Ord Gi* on one numeric attribute over a set of routes (lon/lat + a
# value column). Identical convention to 5b_refugia_hotspot_stats.R's
# run_getis_ord(): k-nearest-neighbor weights with include.self() (true Gi*,
# not plain Gi), BH FDR correction, and the same 7-level confidence class.
run_getis_ord <- function(df, value_col) {
  if (nrow(df) < k_neighbors + 1) {
    message("  Skipping Gi* (need > ", k_neighbors, " routes, have ", nrow(df), ")")
    return(NULL)
  }
  coords <- as.matrix(df[, c("longitude", "latitude")])
  nb <- knearneigh(coords, k = k_neighbors) %>% knn2nb()
  listw <- nb2listw(include.self(nb), style = "B")

  z <- as.numeric(localG(df[[value_col]], listw))
  p <- 2 * pnorm(-abs(z))

  df$gi_z <- z
  df$gi_p <- p
  df$gi_q <- p.adjust(p, method = "BH")
  df$gi_fdr <- dplyr::case_when(
    df$gi_q < fdr_q & z > 0 ~ "Refugia (FDR)",
    df$gi_q < fdr_q & z < 0 ~ "Compound decline (FDR)",
    TRUE                    ~ "Not significant (FDR)"
  )
  df$gi_class <- cut(
    z,
    breaks = c(-Inf, -2.58, -1.96, -1.65, 1.65, 1.96, 2.58, Inf),
    labels = c("Cold Spot (99%)", "Cold Spot (95%)", "Cold Spot (90%)",
              "Not Significant",
              "Hot Spot (90%)", "Hot Spot (95%)", "Hot Spot (99%)")
  )
  df
}

gi_fill_colors <- c(
  "Cold Spot (99%)" = "#2166AC", "Cold Spot (95%)" = "#67A9CF", "Cold Spot (90%)" = "#D1E5F0",
  "Not Significant"  = "grey85",
  "Hot Spot (90%)"   = "#FDDBC7", "Hot Spot (95%)"   = "#EF8A62", "Hot Spot (99%)"  = "#B2182B"
)

make_gi_map <- function(df, title, subtitle, file_path) {
  p <- ggplot() +
    geom_sf(data = us_states, fill = "grey97", color = "grey70", linewidth = 0.3) +
    geom_point(data = df, aes(x = longitude, y = latitude, color = gi_class),
               size = 1.6, alpha = 0.85) +
    geom_point(data = df %>% filter(gi_fdr != "Not significant (FDR)"),
               aes(x = longitude, y = latitude),
               shape = 1, size = 2.8, stroke = 0.5, color = "black") +
    scale_color_manual(values = gi_fill_colors, name = "Hotspot\nclass", drop = FALSE) +
    coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
    labs(title = title, subtitle = subtitle,
         caption = paste0("Red = refugia (locally better despite overall decline + predicted decreasing suitability); ",
                          "blue = compound decline (locally worse). ",
                          "Black ring = still significant after Benjamini-Hochberg FDR correction (q < ", fdr_q, ")"),
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

  # Species' OVERALL trend: mean trend across all of that species' converged
  # routes, any SDM category -- computed once per model_tag (doesn't depend
  # on scenario, since it isn't restricted to a category). Point estimate
  # only, no credible-interval requirement, unless
  # require_credible_for_overall_trend restricts the routes going into this
  # mean to each route's own credible ones first.
  overall_trend_input <- target_sdm %>% filter(!is.na(trend))
  if (require_credible_for_overall_trend) {
    n_before <- nrow(overall_trend_input)
    overall_trend_input <- overall_trend_input %>% filter(trend_lci > 0 | trend_uci < 0)
    cat("require_credible_for_overall_trend = TRUE -> dropped", n_before - nrow(overall_trend_input),
        "of", n_before, "row(s) whose own 90% CI includes zero, before computing species overall trend.\n")
  }

  species_overall <- overall_trend_input %>%
    group_by(species, species_code) %>%
    summarise(overall_trend = mean(trend), n_routes_total = n(), .groups = "drop") %>%
    filter(n_routes_total >= min_routes_overall)

  declining_species <- species_overall %>% filter(overall_trend < 0)
  cat("\n=== model_tag:", model_tag, "===\n")
  cat("Species with overall (all-route) mean trend < 0:", nrow(declining_species),
      "of", nrow(species_overall), "\n")

  write.csv(species_overall,
            file.path(out_dir, paste0(run_label, "_", model_tag,
                                      if (require_credible_for_overall_trend) "_overallcredible" else "",
                                      "_species_overall_trend.csv")),
            row.names = FALSE)

  for (scenario in scenarios) {

    analysis <- target_sdm %>%
      filter(!is.na(.data[[scenario]]), !is.na(trend)) %>%
      transmute(species, species_code, group, route, latitude, longitude, trend,
                trend_lci, trend_uci,
                category = group_category(.data[[scenario]]))

    cat("\n----- model_tag:", model_tag, "| scenario:", scenario, "-----\n")

    # Decreasing-suitability routes for overall-declining species only; the
    # route's own trend sign is NOT filtered -- it stays continuous.
    con <- analysis %>%
      filter(category == "Decreasing suitability", species_code %in% declining_species$species_code)

    if (require_credible_route_trend) {
      n_before <- nrow(con)
      con <- con %>% filter(trend_lci > 0 | trend_uci < 0)
      cat("  require_credible_route_trend = TRUE -> dropped", n_before - nrow(con),
          "of", n_before, "row(s) whose own 90% CI includes zero.\n")
    }

    if (nrow(con) == 0) {
      message("  No qualifying rows -- skipping.")
      next
    }
    cat("  Decreasing-suitability rows for overall-declining species:", nrow(con), "\n")

    route_summary <- con %>%
      group_by(route) %>%
      summarise(latitude = first(latitude), longitude = first(longitude),
                n_species = n(), mean_trend = mean(trend), .groups = "drop")

    cat("  Eligible routes:", nrow(route_summary), "\n")

    gi <- run_getis_ord(route_summary, "mean_trend")
    if (is.null(gi)) next

    n_refugia <- sum(gi$gi_fdr == "Refugia (FDR)")
    n_decline <- sum(gi$gi_fdr == "Compound decline (FDR)")
    cat("  FDR-significant refugia routes (z>0):", n_refugia, "\n")
    cat("  FDR-significant compound-decline routes (z<0):", n_decline, "\n")

    file_prefix <- paste0(run_label, "_", model_tag, "_", scenario,
                          if (length(exclude_species_codes) > 0) "_noinvasive" else "",
                          if (require_credible_for_overall_trend) "_overallcredible" else "",
                          if (require_credible_route_trend) "_credible" else "")

    out_csv <- file.path(out_dir, paste0(file_prefix, "_route_gi.csv"))
    write.csv(gi, out_csv, row.names = FALSE)
    cat("  Wrote:", basename(out_csv), "\n")

    # Per-route species detail behind the FDR-significant routes in each
    # direction, same shape as 5b's fdr_species outputs.
    refugia_species <- con %>% filter(route %in% gi$route[gi$gi_fdr == "Refugia (FDR)"]) %>%
      count(route, species, group, sort = TRUE)
    write.csv(refugia_species, file.path(out_dir, paste0(file_prefix, "_refugia_species.csv")),
              row.names = FALSE)

    decline_species <- con %>% filter(route %in% gi$route[gi$gi_fdr == "Compound decline (FDR)"]) %>%
      count(route, species, group, sort = TRUE)
    write.csv(decline_species, file.path(out_dir, paste0(file_prefix, "_compound_decline_species.csv")),
              row.names = FALSE)

    model_scenario_txt <- paste0(model_tag, ", ", toupper(sub("rcp", "RCP ", scenario)))
    make_gi_map(
      gi,
      paste0("Overall-declining species, Decreasing-suitability routes (", model_scenario_txt, ")"),
      paste0(nrow(gi), " eligible routes | ", n_refugia, " refugia (FDR) | ", n_decline, " compound decline (FDR)"),
      file.path(plot_dir, paste0(file_prefix, "_gi_map.png"))
    )
  }
}

cat("\n=== Overall-decline hotspot statistics complete (model_tags = ", paste(model_tags, collapse = ", "),
    ", scenarios = ", paste(scenarios, collapse = ", "), ", bird_group = ",
    if (is.na(bird_group)) "all_groups (NA)" else bird_group, ") ===\n", sep = "")
