# 5b_refugia_hotspot_stats.R
#
# Follow-up to 5_hotspots.R. 5_hotspots.R's pooled, all-species mismatch RATE
# maps (grid heatmap / bubble map) look diffuse -- no obvious regional
# cluster -- which is exactly what you'd expect if either (a) there's no real
# spatial clustering and the mismatches are just scattered noise, or (b) real
# local clusters exist but are being averaged away by eye when looking at a
# rate map, or by pooling ecologically dissimilar species together. This
# script asks the same underlying question -- where is climate-suitability
# CONTRACTION met with a species actually persisting/increasing, i.e.
# candidate refugia -- three more targeted ways:
#
#   PART A -- Getis-Ord Gi* hotspot statistic on route-level mismatch RATE.
#     A formal local-spatial-autocorrelation test (via spdep), not a visual
#     read of a heatmap: for each route, is its mismatch rate significantly
#     higher (or lower) than chance given its k-nearest neighbors' rates?
#     This can surface statistically real local clusters that a smoothed or
#     binned map hides inside a lot of scattered noise, or conversely can
#     confirm that "it really is diffuse" is the correct read.
#
#   PART B -- the same Gi* statistic, but on a MAGNITUDE-weighted "refugia
#     strength index" (mean trend among a route's credibly-mismatching
#     species) instead of a 0/1 rate. A route where mismatching species blow
#     past their SDM prediction by +8%/year is a stronger refugia candidate
#     than one that just barely clears trend_lci > 0 -- the rate-only view
#     in Part A (and in 5_hotspots.R) can't distinguish those. Computed only
#     on routes with >= 1 mismatching species (undefined otherwise), with
#     spatial neighbors taken among that same restricted set.
#
#   PART C -- cross-species/cross-guild convergence. A route where several
#     species from DIFFERENT habitat groups (data/spp_names_codes_group_aou.csv's
#     Group, carried through as "group" here) independently buck their own
#     SDM prediction is stronger refugia evidence than a route where only one
#     species (or several from the same guild, who may share one idiosyncratic
#     reason) does -- convergent evidence across independent lines beats a
#     single line of evidence. Reports n_groups_mismatch per route alongside
#     Part A/B's statistics rather than as its own spatial model.
#
# Multiple testing: thousands of routes are tested at once, so the per-test
# 90/95/99% classes (gi_class_*) overstate how many hot spots are real. Every
# analysis also gets a Benjamini-Hochberg false-discovery-rate adjusted
# q-value (gi_q_*) and a gi_fdr_* flag (q < fdr_q) -- the FDR-surviving
# routes are the ones to report. 5c_gi_k_sensitivity.R adds robustness to the
# neighborhood-size choice on top of this.
#
# Only the "credible" mismatch definition (trend_lci > 0 / trend_uci < 0,
# see 5_hotspots.R) is used here -- the project settled on credibility over a
# fixed magnitude threshold as the primary criterion; see the project
# record/methods note for the full reasoning.
#
# NOT attempted here (documented, not implemented): cross-referencing
# hotspot routes against the raw `anthro` land-use covariate (data/Anthro.csv)
# to test whether refugia coincide with low anthropogenic pressure -- that
# raw per-route-year covariate file isn't part of this pipeline's committed
# output (only the fitted models' summaries are), so this would need
# data/Anthro.csv to be available locally. Worth doing if/when it is.
#
# Reads:  output/species_routes_covariates/per_species_sdm/*_route_trends_sdm.csv
#         (written by 3c_add_SDM_covariates.R) -- same input as 5_hotspots.R
# Writes: output/species_routes_covariates/hotspots/refugia_stats/<run_label>_<model_tag>_<scenario>_<section>_route_gi.csv
#         output/species_routes_covariates/hotspots/refugia_stats/plots/<run_label>_<model_tag>_<scenario>_<section>_gi_rate_map.png
#         output/species_routes_covariates/hotspots/refugia_stats/plots/<run_label>_<model_tag>_<scenario>_<section>_gi_magnitude_map.png

library(here)
library(tidyverse)
library(sf)
library(spdep)

here::i_am("5b_refugia_hotspot_stats.R")

# Settings -------------------------------------------------------------------
bird_group <- NA   # NA pools every group (matches 5_hotspots.R's default);
                   # set to one Group value to restrict, e.g. to test whether
                   # pooling guilds together is washing out a per-guild signal.

model_tags <- c("base", "anthro")
scenarios  <- c("rcp45", "rcp85")

require_route_converged <- TRUE

# A route needs at least this many contraction/expansion-predicted species
# before its mismatch RATE is included in Part A's Gi* -- otherwise a route
# with n_category = 1 contributes a meaningless 0%/100% extreme into the
# spatial statistic. Same rationale as 5_hotspots.R's min_rows_per_cell, just
# at the route level instead of the grid-cell level.
min_n_category <- 3

# Number of nearest-neighbor routes used to build the spatial weights matrix
# for Getis-Ord Gi*. Routes aren't on a regular grid, so k-nearest-neighbors
# (not a fixed distance band) is used -- keeps each route's neighborhood size
# comparable even where route density varies a lot (dense in the East,
# sparse in the interior West).
k_neighbors <- 8

# Benjamini-Hochberg FDR threshold: of the routes flagged as FDR-significant,
# expect at most this share to be false positives.
fdr_q <- 0.05

in_dir  <- here::here("output", "species_routes_covariates", "per_species_sdm")
out_dir <- here::here("output", "species_routes_covariates", "hotspots", "refugia_stats")
plot_dir <- file.path(out_dir, "plots")
if (!dir.exists(out_dir))  dir.create(out_dir,  recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Read every SDM CSV (same input as 5_hotspots.R) -----------------------------
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

if (nrow(all_sdm_raw) == 0) {
  stop("No rows left after filtering (bird_group = '", bird_group,
       "', require_route_converged = ", require_route_converged, ").")
}

# Contraction = 1,2,3 | Stable = 4 | Expansion = 5,6,7 -- kept in sync with
# 4c_statistical_analysis_and_visualization_covariates.R / 5_hotspots.R.
group_category <- function(x) {
  dplyr::case_when(
    x %in% c(1, 2, 3) ~ "Contraction",
    x == 4            ~ "Stable",
    x %in% c(5, 6, 7) ~ "Expansion",
    TRUE              ~ NA_character_
  )
}

# "credible" mismatch sections only -- see header comment for why.
mismatch_definitions <- list(
  list(key = "contraction_credibly_increasing",
       category = "Contraction",
       test = function(d) d$trend_lci > 0,
       label = "Contraction-predicted but credibly increasing")
  # Expansion section paused to focus on contraction; restore by adding a
  # comma after the entry above and uncommenting:
  # list(key = "expansion_credibly_decreasing",
  #      category = "Expansion",
  #      test = function(d) d$trend_uci < 0,
  #      label = "Expansion-predicted but credibly decreasing")
)

# US-states background map (same as 5_hotspots.R).
data("us_states", package = "spData")
us_states <- st_transform(us_states, 4326)
us_bbox <- st_bbox(us_states)
map_xlim <- c(us_bbox["xmin"] - 1, us_bbox["xmax"] + 1)
map_ylim <- c(us_bbox["ymin"] - 1, us_bbox["ymax"] + 1)

# Getis-Ord Gi* on one numeric attribute over a set of routes (lon/lat + a
# value column). Builds k-nearest-neighbor spatial weights, runs localG, and
# classifies each route into the conventional 7-level confidence bin (same
# convention as ArcGIS's Hot Spot Analysis tool: |z| >= 2.58 -> 99%,
# >= 1.96 -> 95%, >= 1.65 -> 90%, else not significant).
run_getis_ord <- function(df, value_col) {
  if (nrow(df) < k_neighbors + 1) {
    message("  Skipping Gi* (need > ", k_neighbors, " routes, have ", nrow(df), ")")
    return(NULL)
  }
  coords <- as.matrix(df[, c("longitude", "latitude")])
  nb <- knearneigh(coords, k = k_neighbors) %>% knn2nb()
  # include.self() counts each route in its own neighborhood, which is what
  # makes spdep::localG return Gi* (without it, it returns plain Gi).
  listw <- nb2listw(include.self(nb), style = "B")   # binary weights

  z <- as.numeric(localG(df[[value_col]], listw))
  p <- 2 * pnorm(-abs(z))   # two-sided p-value from the standard normal z

  df$gi_z <- z
  df$gi_p <- p
  # Benjamini-Hochberg false-discovery-rate adjustment across every route
  # tested in THIS analysis (the per-test 90/95/99% classes below don't
  # account for running thousands of tests at once). gi_fdr flags routes
  # whose adjusted q < fdr_q, split by direction.
  df$gi_q <- p.adjust(p, method = "BH")
  df$gi_fdr <- dplyr::case_when(
    df$gi_q < fdr_q & z > 0 ~ "Hot (FDR)",
    df$gi_q < fdr_q & z < 0 ~ "Cold (FDR)",
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

make_gi_map <- function(df, fill_var, title, subtitle, legend_name, file_path) {
  p <- ggplot() +
    geom_sf(data = us_states, fill = "grey97", color = "grey70", linewidth = 0.3) +
    geom_point(data = df, aes(x = longitude, y = latitude, color = .data[[fill_var]]),
               size = 1.6, alpha = 0.85) +
    # Black ring = also significant after false-discovery-rate correction.
    geom_point(data = df %>% filter(gi_fdr != "Not significant (FDR)"),
               aes(x = longitude, y = latitude),
               shape = 1, size = 2.8, stroke = 0.5, color = "black") +
    scale_color_manual(values = gi_fill_colors, name = legend_name, drop = FALSE) +
    coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
    labs(title = title, subtitle = subtitle,
         caption = paste0("Black ring = still significant after Benjamini-Hochberg FDR correction (q < ",
                          fdr_q, ")"),
         x = NULL, y = NULL) +
    theme_minimal() +
    theme(plot.title    = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 10),
          plot.caption  = element_text(size = 8, color = "grey40"),
          axis.text     = element_text(size = 8))
  ggsave(file_path, p, width = 9, height = 6.5, dpi = 150)
  cat("Saved plot:", basename(file_path), "\n")
}

# Main loop: model_tag x scenario x mismatch-section --------------------------
for (model_tag in model_tags) {
  target_sdm <- all_sdm_raw %>% filter(model == model_tag)
  if (nrow(target_sdm) == 0) {
    message("Skipping model_tag = '", model_tag, "' -- no rows found.")
    next
  }

  for (scenario in scenarios) {

    analysis <- target_sdm %>%
      filter(!is.na(.data[[scenario]]), !is.na(trend), !is.na(trend_lci), !is.na(trend_uci)) %>%
      transmute(species, species_code, group, route, latitude, longitude,
                trend, trend_lci, trend_uci,
                category = group_category(.data[[scenario]]))

    for (def in mismatch_definitions) {
      cat("\n----- model_tag:", model_tag, "| scenario:", scenario,
          "| section:", def$key, "-----\n")

      category_df <- analysis %>% filter(category == def$category)
      category_df$mismatch <- def$test(category_df)

      if (nrow(category_df) == 0) {
        message("  No ", def$category, "-category rows -- skipping ", def$key, ".")
        next
      }

      # Route-level summary, with Part B/C additions beyond 5_hotspots.R's
      # route_summary: mean_trend_mismatch (Part B's magnitude index) and
      # n_groups_mismatch (Part C's cross-guild convergence count).
      route_summary <- category_df %>%
        group_by(route) %>%
        summarise(
          latitude  = first(latitude), longitude = first(longitude),
          n_category = n(), n_mismatch = sum(mismatch),
          pct_mismatch = n_mismatch / n_category,
          mean_trend_mismatch = if (any(mismatch)) mean(trend[mismatch]) else NA_real_,
          n_groups_mismatch = n_distinct(group[mismatch]),
          .groups = "drop"
        )

      cat("  Routes:", nrow(route_summary),
          "| with >=1 mismatch:", sum(route_summary$n_mismatch > 0), "\n")

      file_prefix <- paste0(run_label, "_", model_tag, "_", scenario, "_", def$key)
      model_scenario_txt <- paste0(model_tag, ", ", toupper(sub("rcp", "RCP ", scenario)))
      subtitle_base <- paste0(def$label, " (", model_scenario_txt, ")")

      # PART A -- Gi* on mismatch RATE, routes with enough sample to trust
      # a rate at all (min_n_category).
      rate_input <- route_summary %>% filter(n_category >= min_n_category)
      cat("  Part A (rate Gi*): ", nrow(rate_input), " routes meet n_category >= ",
          min_n_category, "\n", sep = "")
      rate_gi <- run_getis_ord(rate_input, "pct_mismatch")

      # PART B -- Gi* on magnitude (mean trend among mismatching species),
      # restricted to routes where that's even defined (>=1 mismatch), with
      # neighbors taken among that same restricted set (see header comment).
      magnitude_input <- route_summary %>% filter(n_mismatch >= 1)
      cat("  Part B (magnitude Gi*): ", nrow(magnitude_input),
          " routes have >=1 mismatching species\n", sep = "")
      magnitude_gi <- run_getis_ord(magnitude_input, "mean_trend_mismatch")

      # Write combined route-level output (Part A + B + C columns together;
      # gi_z/gi_p/gi_class get an "_rate"/"_magnitude" suffix to keep the two
      # Gi* runs distinct in one file).
      out_df <- route_summary
      if (!is.null(rate_gi)) {
        out_df <- out_df %>%
          left_join(rate_gi %>% select(route, gi_z_rate = gi_z, gi_p_rate = gi_p,
                                       gi_q_rate = gi_q, gi_class_rate = gi_class,
                                       gi_fdr_rate = gi_fdr),
                    by = "route")
      }
      if (!is.null(magnitude_gi)) {
        out_df <- out_df %>%
          left_join(magnitude_gi %>% select(route, gi_z_magnitude = gi_z, gi_p_magnitude = gi_p,
                                            gi_q_magnitude = gi_q, gi_class_magnitude = gi_class,
                                            gi_fdr_magnitude = gi_fdr),
                    by = "route")
      }
      out_csv <- file.path(out_dir, paste0(file_prefix, "_route_gi.csv"))
      write.csv(out_df, out_csv, row.names = FALSE)
      cat("  Wrote:", basename(out_csv), "\n")

      # Part C summary (printed, not a separate file): does higher mismatch
      # rate coexist with more cross-guild convergence, or is a high rate
      # often just one guild repeating itself?
      conv_summary <- route_summary %>% filter(n_mismatch >= 1) %>%
        summarise(n_routes = n(),
                  pct_multi_guild = 100 * mean(n_groups_mismatch >= 2),
                  mean_n_groups = mean(n_groups_mismatch))
      cat("  Part C: of routes with >=1 mismatch, ",
          sprintf("%.1f%%", conv_summary$pct_multi_guild),
          " involve >=2 different habitat groups (mean groups/route = ",
          sprintf("%.2f", conv_summary$mean_n_groups), ")\n", sep = "")

      # Plots -------------------------------------------------------------
      if (!is.null(rate_gi)) {
        n_hot_rate <- sum(grepl("^Hot", rate_gi$gi_class))
        n_cold_rate <- sum(grepl("^Cold", rate_gi$gi_class))
        n_fdr_rate <- sum(rate_gi$gi_fdr == "Hot (FDR)")
        make_gi_map(
          rate_gi, "gi_class",
          paste0("Gi* hotspots: mismatch rate (", model_scenario_txt, ")"),
          paste0(subtitle_base, "\n", nrow(rate_gi), " routes (n_category >= ", min_n_category, ") | ",
                n_hot_rate, " hot-spot (", n_fdr_rate, " after FDR), ", n_cold_rate,
                " cold-spot routes"),
          "Rate hotspot\nclass",
          file.path(plot_dir, paste0(file_prefix, "_gi_rate_map.png"))
        )
      }

      if (!is.null(magnitude_gi)) {
        n_hot_mag <- sum(grepl("^Hot", magnitude_gi$gi_class))
        n_cold_mag <- sum(grepl("^Cold", magnitude_gi$gi_class))
        n_fdr_mag <- sum(magnitude_gi$gi_fdr == "Hot (FDR)")
        make_gi_map(
          magnitude_gi, "gi_class",
          paste0("Gi* hotspots: refugia strength (", model_scenario_txt, ")"),
          paste0(subtitle_base, "\n", nrow(magnitude_gi), " routes with >=1 mismatch | ",
                n_hot_mag, " hot-spot (", n_fdr_mag, " after FDR), ", n_cold_mag,
                " cold-spot routes"),
          "Magnitude\nhotspot class",
          file.path(plot_dir, paste0(file_prefix, "_gi_magnitude_map.png"))
        )
      }
    }
  }
}

cat("\n=== Refugia hotspot statistics complete (model_tags = ", paste(model_tags, collapse = ", "),
    ", scenarios = ", paste(scenarios, collapse = ", "), ", bird_group = ",
    if (is.na(bird_group)) "all_groups (NA)" else bird_group, ") ===\n", sep = "")
