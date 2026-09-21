# 5_hotspots.R
#
# Maps route x species instances where the SDM classified-change model and
# the route-level population trend from the iCAR model DISAGREE -- two such
# mismatches ("sections" below), each run for every model_tag x scenario
# combination. Two mismatch-definition MODES are supported (mismatch_mode
# setting below), each writing to its own subfolder so runs don't collide:
#
#   "credible"  -- a route only counts if its 90% posterior credible
#                  interval on trend (trend_lci/trend_uci, from 2c) excludes
#                  zero in the direction that disagrees with the SDM. This
#                  is the same "credible effect" convention the pipeline
#                  already uses for gamma1 (gamma1_excludes_zero in 2c),
#                  applied here to trend instead.
#   "threshold" -- a route counts if its point-estimate trend magnitude
#                  clears a fixed %/year cutoff (threshold_magnitude below,
#                  mirrored for both directions), regardless of how wide its
#                  credible interval is.
#
# See the project record/methods note for the full reasoning: a fixed
# %/year magnitude threshold answers "is this ecologically meaningful,"
# while credibility answers "do we trust this specific estimate" -- two
# different questions the model already lets us ask separately, which is
# why both modes are kept available here rather than picking one.
#
#   1. contraction_credibly_increasing / contraction_increasing_<X>pct --
#      SDM predicts CONTRACTION but the route trend is credibly increasing
#      (trend_lci > 0) / increasing by at least X%/year (trend >= X).
#   2. expansion_credibly_decreasing / expansion_decreasing_<X>pct -- the
#      mirror-image mismatch: SDM predicts EXPANSION but the route trend is
#      credibly decreasing (trend_uci < 0) / decreasing by at least X%/year
#      (trend < -X).
#
# Both sections ask the same underlying question -- where does the climate-
# suitability model's directional prediction disagree with what the route
# is actually doing -- just for opposite (predicted direction, observed
# direction) pairs. These are candidate "hotspots": places worth a closer
# look (lagged response, a covariate the SDM misses, local refugia, etc.).
#
# Contraction / Stable / Expansion category boundaries are the SAME grouping
# 4c_statistical_analysis_and_visualization_covariates.R uses (Contraction =
# rcp category 1,2,3 | Stable = 4 | Expansion = 5,6,7); keep group_category()
# below in sync with that script's copy if the boundaries ever change.
#
# Three plots per model_tag x scenario x mismatch-section combination, all
# on a US-states map background (contiguous US only -- spData::us_states --
# which matches this project's route coverage):
#   1. Scatter map: one point per mismatching (species, route) pair, colored
#      by route trend. Heavy overlap of points at one location is itself a
#      visual signal of a hotspot, but this is raw point density -- an area
#      with more BBS routes shows up denser even if its mismatch RATE isn't
#      actually higher.
#   2. Grid heatmap: routes binned into a fixed-size lon/lat grid, each cell
#      colored by its own pct_mismatch (mismatches / SDM-category-predicted
#      species pooled across every route in that cell). Rate-based and local
#      by construction -- no cross-region smoothing, unlike a KDE surface --
#      but the bin size (grid_bin_size_deg) is an explicit resolution choice,
#      and cells with too few predicted-category rows are dropped
#      (min_rows_per_cell) so a single noisy route can't paint a whole cell
#      red/yellow.
#   3. Route bubble map: one point per ROUTE (not per species-route pair,
#      unlike plot 1) at its exact coordinates, colored by that route's own
#      pct_mismatch and sized by n_category (how many SDM-category-predicted
#      species were evaluated there, i.e. how reliable that route's rate
#      is). The finest-grained, least-smoothed view -- no binning or density
#      estimation at all.
#
# Reads:  output/species_routes_covariates/per_species_sdm/*_route_trends_sdm.csv
#         (written by 3c_add_SDM_covariates.R)
# Writes (per model_tag x scenario x mismatch-section), under a subfolder
# named after mismatch_mode ("credible" or "threshold") so the two modes'
# outputs never collide:
#         output/species_routes_covariates/hotspots/<mismatch_mode>/<run_label>_<model_tag>_<scenario>_<section>_mismatch_routes.csv
#         output/species_routes_covariates/hotspots/<mismatch_mode>/<run_label>_<model_tag>_<scenario>_<section>_route_summary.csv
#         output/species_routes_covariates/hotspots/<mismatch_mode>/plots/<run_label>_<model_tag>_<scenario>_<section>_mismatch_scatter_map.png
#         output/species_routes_covariates/hotspots/<mismatch_mode>/plots/<run_label>_<model_tag>_<scenario>_<section>_mismatch_grid_heatmap.png
#         output/species_routes_covariates/hotspots/<mismatch_mode>/plots/<run_label>_<model_tag>_<scenario>_<section>_mismatch_bubble_map.png

library(here)
library(tidyverse)
library(sf)

here::i_am("5_hotspots.R")

# Settings -------------------------------------------------------------------
bird_group <- NA   # one of the 12 Group values in
                   # data/spp_names_codes_group_aou.csv, or NA to pool every
                   # group together (matches 4c's default).

model_tags <- c("base", "anthro")   # loop over both, as 4c does.
scenarios  <- c("rcp45", "rcp85")   # loop over both climate scenarios.

# "credible" -- mismatch = 90% CI excludes zero (trend_lci/trend_uci).
# "threshold" -- mismatch = point-estimate trend magnitude >= threshold_magnitude
# (mirrored for both directions). See the header comment for the reasoning
# behind keeping both modes available rather than picking one.
mismatch_mode <- "credible"

# Only used when mismatch_mode == "threshold": %/year. Applied as
# trend >= +threshold_magnitude (contraction section, "increasing") and
# trend < -threshold_magnitude (expansion section, "decreasing"). 1.2 is
# PIF's Watch List "50% cumulative loss since the mid-1960s BBS baseline"
# criterion, annualized (1 - 0.5^(1/58) ~= 1.2%/year) and mirrored to the
# increasing direction rather than combining it with a differently-derived
# benchmark for the other direction -- see the project record/methods note.
threshold_magnitude <- 1.2

# Same convention as 4c: drop routes whose own alpha[r]/beta[r] didn't
# individually meet Rhat < 1.01 & bulk ESS > 400 before mapping anything.
require_route_converged <- TRUE

# Grid heatmap resolution (degrees lon/lat per cell) and the minimum number
# of SDM-category-predicted (species, route) rows a cell needs before its
# pct_mismatch is shown at all -- below this it's dropped (grey) rather than
# plotted, since e.g. 1/1 = 100% from a single route is not a rate estimate.
grid_bin_size_deg   <- 1.5
min_rows_per_cell   <- 5

if (!mismatch_mode %in% c("credible", "threshold")) {
  stop("mismatch_mode must be 'credible' or 'threshold', got: '", mismatch_mode, "'")
}

in_dir  <- here::here("output", "species_routes_covariates", "per_species_sdm")
out_dir <- here::here("output", "species_routes_covariates", "hotspots", mismatch_mode)
plot_dir <- file.path(out_dir, "plots")
if (!dir.exists(out_dir))  dir.create(out_dir,  recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Read every SDM CSV (same input as 4c) ---------------------------------------
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

# Contraction = 1,2,3 | Stable = 4 | Expansion = 5,6,7 | anything else
# (0 = never suitable, or an out-of-legend raster value) -> NA, excluded.
# Kept in sync with 4c_statistical_analysis_and_visualization_covariates.R.
group_category <- function(x) {
  dplyr::case_when(
    x %in% c(1, 2, 3) ~ "Contraction",
    x == 4            ~ "Stable",
    x %in% c(5, 6, 7) ~ "Expansion",
    TRUE              ~ NA_character_
  )
}

# The two mismatch sections described in the header comment, built for
# whichever mismatch_mode is active. Each ties an SDM category to a test
# (on the filtered, category-matched tibble, so it can reach trend_lci/
# trend_uci as well as trend) and to the labels used in filenames/titles/
# legends below -- adding another section only requires one more entry,
# the loop below is generic.
if (mismatch_mode == "credible") {
  mismatch_definitions <- list(
    list(key = "contraction_credibly_increasing",
         category = "Contraction",
         test = function(d) d$trend_lci > 0,
         label = "Contraction-predicted but credibly increasing",
         rate_label = "Credibly increasing rate\n(of contraction-\npredicted species)")
    # Expansion section paused to focus on contraction; restore by adding a
    # comma after the entry above and uncommenting:
    # list(key = "expansion_credibly_decreasing",
    #      category = "Expansion",
    #      test = function(d) d$trend_uci < 0,
    #      label = "Expansion-predicted but credibly decreasing",
    #      rate_label = "Credibly decreasing rate\n(of expansion-\npredicted species)")
  )
} else {
  mismatch_definitions <- list(
    list(key = paste0("contraction_increasing_", threshold_magnitude, "pct"),
         category = "Contraction",
         test = function(d) d$trend >= threshold_magnitude,
         label = paste0("Contraction-predicted but increasing ≥", threshold_magnitude, "%/yr"),
         rate_label = paste0("≥", threshold_magnitude, "%/yr increasing rate\n(of contraction-\npredicted species)"))
    # Expansion section paused to focus on contraction; restore by adding a
    # comma after the entry above and uncommenting:
    # list(key = paste0("expansion_decreasing_", threshold_magnitude, "pct"),
    #      category = "Expansion",
    #      test = function(d) d$trend < -threshold_magnitude,
    #      label = paste0("Expansion-predicted but decreasing ≥", threshold_magnitude, "%/yr"),
    #      rate_label = paste0("≥", threshold_magnitude, "%/yr decreasing rate\n(of expansion-\npredicted species)"))
  )
}

# US-states background map, contiguous US only (matches this project's route
# coverage: lower 48, no Alaska/Hawaii/Canada) -- bundled with the already-
# installed spData package, so no download/API call needed.
data("us_states", package = "spData")
us_states <- st_transform(us_states, 4326)
us_bbox <- st_bbox(us_states)
map_xlim <- c(us_bbox["xmin"] - 1, us_bbox["xmax"] + 1)
map_ylim <- c(us_bbox["ymin"] - 1, us_bbox["ymax"] + 1)

# Shown as a small caption on every plot rather than folded into the
# subtitle, which is already busy with route/species counts and would
# otherwise run past the plot edge.
mismatch_caption <- if (mismatch_mode == "credible") {
  "Mismatch criterion: route's 90% posterior credible interval on trend excludes zero"
} else {
  paste0("Mismatch criterion: point-estimate trend magnitude ≥ ", threshold_magnitude,
        "%/year (fixed threshold; credible interval not required)")
}

# Scatter map: one point per mismatching (species, route) pair -------------
make_scatter_map <- function(pts, title, subtitle, file_path) {
  # Clip both tails at the 5th/95th percentile rather than assuming trend is
  # one-signed -- section 1 mismatches are bounded below (>= 0-ish) with
  # outliers on the high side, section 2 (decreasing) is bounded above (< 0)
  # with outliers on the low side. Percentile clipping on both ends handles
  # either shape without special-casing per section.
  color_lims <- quantile(pts$trend, c(0.05, 0.95), na.rm = TRUE)
  p <- ggplot() +
    geom_sf(data = us_states, fill = "grey96", color = "grey65", linewidth = 0.3) +
    geom_point(data = pts, aes(x = longitude, y = latitude, color = trend),
               alpha = 0.45, size = 1.8) +
    scale_color_viridis_c(option = "C", name = "Route trend\n(%/yr)",
                          limits = color_lims, oob = scales::squish) +
    coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
    labs(title = title, subtitle = subtitle, caption = mismatch_caption, x = NULL, y = NULL) +
    theme_minimal() +
    theme(plot.title    = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 11),
          plot.caption  = element_text(size = 8, color = "grey40"),
          axis.text     = element_text(size = 8))
  ggsave(file_path, p, width = 9, height = 6.5, dpi = 150)
  cat("Saved plot:", basename(file_path), "\n")
}

# Grid heatmap: bin ALL SDM-category-predicted rows (mismatch TRUE or FALSE)
# into fixed-size lon/lat cells and color each cell by its own pct_mismatch.
# Local and rate-based -- the opposite failure mode from a KDE surface,
# which shows raw point density (confounded with route density) smoothed
# across whatever bandwidth you pick.
make_grid_heatmap <- function(category_pts, rate_label, title, subtitle, file_path) {
  cell_summary <- category_pts %>%
    mutate(lon_bin = floor(longitude / grid_bin_size_deg) * grid_bin_size_deg,
           lat_bin = floor(latitude  / grid_bin_size_deg) * grid_bin_size_deg) %>%
    group_by(lon_bin, lat_bin) %>%
    summarise(n_category = n(), n_mismatch = sum(mismatch),
              pct_mismatch = n_mismatch / n_category, .groups = "drop") %>%
    filter(n_category >= min_rows_per_cell)

  if (nrow(cell_summary) == 0) {
    message("  Skipping grid heatmap (no cell reached min_rows_per_cell = ",
            min_rows_per_cell, "): ", basename(file_path))
    return(invisible(NULL))
  }

  p <- ggplot() +
    geom_sf(data = us_states, fill = "grey97", color = "grey70", linewidth = 0.3) +
    geom_tile(data = cell_summary,
              aes(x = lon_bin + grid_bin_size_deg / 2, y = lat_bin + grid_bin_size_deg / 2,
                  fill = pct_mismatch),
              width = grid_bin_size_deg, height = grid_bin_size_deg, alpha = 0.85) +
    scale_fill_viridis_c(option = "C", name = rate_label,
                         labels = scales::percent, limits = c(0, 1)) +
    geom_sf(data = us_states, fill = NA, color = "grey30", linewidth = 0.3) +
    coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
    labs(title = title, subtitle = subtitle,
         caption = paste0(mismatch_caption, "\n", grid_bin_size_deg, "° cells, n ≥ ",
                          min_rows_per_cell, " predicted-category rows each"),
         x = NULL, y = NULL) +
    theme_minimal() +
    theme(plot.title    = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 10),
          plot.caption  = element_text(size = 8, color = "grey40"),
          axis.text     = element_text(size = 8))
  ggsave(file_path, p, width = 9, height = 6.5, dpi = 150)
  cat("Saved plot:", basename(file_path), "\n")
}

# Route bubble map: one point per route (route_summary, already aggregated),
# colored by that route's own pct_mismatch and sized by n_category (how many
# SDM-category-predicted species were evaluated there) as a visual cue for
# how much to trust that route's rate. No binning or smoothing at all -- the
# most local view of the three.
make_bubble_map <- function(route_summary, rate_label, title, subtitle, file_path) {
  p <- ggplot() +
    geom_sf(data = us_states, fill = "grey96", color = "grey65", linewidth = 0.3) +
    geom_point(data = route_summary,
               aes(x = longitude, y = latitude, color = pct_mismatch, size = n_category),
               alpha = 0.75) +
    scale_color_viridis_c(option = "C", name = rate_label,
                          labels = scales::percent, limits = c(0, 1)) +
    scale_size_continuous(name = "Species\nevaluated", range = c(0.6, 4)) +
    coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
    labs(title = title, subtitle = subtitle, caption = mismatch_caption, x = NULL, y = NULL) +
    theme_minimal() +
    theme(plot.title    = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 11),
          plot.caption  = element_text(size = 8, color = "grey40"),
          axis.text     = element_text(size = 8))
  ggsave(file_path, p, width = 9, height = 6.5, dpi = 150)
  cat("Saved plot:", basename(file_path), "\n")
}

# Main loop: model_tag x scenario x mismatch-section -------------------------
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
        message("  No ", def$category, "-category rows for ", model_tag, "/", scenario,
                " -- skipping ", def$key, ".")
        next
      }

      mismatch_rows <- category_df %>% filter(mismatch)

      cat("  ", def$category, "-predicted (species, route) rows: ", nrow(category_df), "\n", sep = "")
      cat("  Of those, mismatching (", def$key, ") rows: ", nrow(mismatch_rows),
          sprintf(" (%.1f%%)\n", 100 * nrow(mismatch_rows) / nrow(category_df)), sep = "")

      file_prefix <- paste0(run_label, "_", model_tag, "_", scenario, "_", def$key)

      mismatch_csv <- file.path(out_dir, paste0(file_prefix, "_mismatch_routes.csv"))
      write.csv(mismatch_rows, mismatch_csv, row.names = FALSE)
      cat("  Wrote:", basename(mismatch_csv), "\n")

      route_summary <- category_df %>%
        group_by(route) %>%
        summarise(latitude = first(latitude), longitude = first(longitude),
                  n_category = n(), n_mismatch = sum(mismatch),
                  pct_mismatch = n_mismatch / n_category, .groups = "drop") %>%
        arrange(desc(pct_mismatch), desc(n_mismatch))

      route_summary_csv <- file.path(out_dir, paste0(file_prefix, "_route_summary.csv"))
      write.csv(route_summary, route_summary_csv, row.names = FALSE)
      cat("  Wrote:", basename(route_summary_csv), "\n")

      if (nrow(mismatch_rows) == 0) {
        message("  No mismatch (species, route) rows for ", model_tag, "/", scenario,
                "/", def$key, " -- skipping plots.")
        next
      }

      title_txt <- paste0(def$label, " routes (",
                          model_tag, ", ", toupper(sub("rcp", "RCP ", scenario)), ")")
      subtitle_txt <- paste0(nrow(mismatch_rows), " (species, route) mismatches of ",
                             nrow(category_df), " ", tolower(def$category), "-predicted, across ",
                             length(unique(mismatch_rows$route)), " routes and ",
                             length(unique(mismatch_rows$species_code)), " species (",
                             run_label, ")")

      make_scatter_map(
        mismatch_rows, title_txt, subtitle_txt,
        file.path(plot_dir, paste0(file_prefix, "_mismatch_scatter_map.png"))
      )

      make_grid_heatmap(
        category_df, def$rate_label, title_txt, subtitle_txt,
        file.path(plot_dir, paste0(file_prefix, "_mismatch_grid_heatmap.png"))
      )

      make_bubble_map(
        route_summary, def$rate_label, title_txt, subtitle_txt,
        file.path(plot_dir, paste0(file_prefix, "_mismatch_bubble_map.png"))
      )
    }
  }
}

cat("\n=== Hotspot mapping complete (model_tags = ", paste(model_tags, collapse = ", "),
    ", scenarios = ", paste(scenarios, collapse = ", "), ", sections = ",
    paste(map_chr(mismatch_definitions, "key"), collapse = ", "), ", bird_group = ",
    if (is.na(bird_group)) "all_groups (NA)" else bird_group, ") ===\n", sep = "")
