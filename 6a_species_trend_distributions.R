# 6a_species_trend_distributions.R
#
# One "raincloud" figure per species x model (base / anthro) x climate
# scenario (RCP 4.5 / RCP 8.5): for every route belonging to that species,
# pulls its population trend and its SDM range-shift category (Contraction /
# Stable / Expansion, from the same classified-change raster values 3d's
# maps use), then shows:
#   - a density curve of route trend PER category (Contraction / Stable /
#     Expansion), filled with the same area colors as
#     3d_visualization_map.R's range-shift maps, so this figure and that map
#     read as one system;
#   - a dashed grey OUTLINE-ONLY density behind all three, for every
#     categorized route regardless of category (the species' overall trend
#     distribution, deliberately drawn with no fill so it never competes
#     with the three colored areas);
#   - a jittered strip of individual route trends below the curves, one row
#     per category, colored the same way;
#   - a black open ring on any route whose OWN 90% credible interval
#     (trend_lci/trend_uci) excludes zero -- same "still significant" ring
#     convention 5b_refugia_hotspot_stats.R uses for FDR-significant
#     hotspots, applied here to route-level credibility instead;
#   - a filled triangle at each category's mean trend.
#
# Companion to 3d_visualization_map.R (same category colors, same route
# trend/credible-interval source) and 4c_statistical_analysis_and_visualization_covariates.R
# (same Contraction/Stable/Expansion grouping and per-species scope) -- this
# figure is the "pull every route's trend out of the map and look at its
# distribution" complement to 3d's spatial view.
#
# Reads:  output/species_routes_covariates/per_species_sdm/*_route_trends_sdm.csv
#           (3c; species, group, model, route, trend, trend_lci/uci,
#           route_converged, rcp45, rcp85 -- rcp45/rcp85 are the raw 0-7
#           classified-change raster codes at that route's location)
# Writes: output/species_trend_distributions/<model>/<scenario>/<species>_<model>_<scenario>_<firstYear>_<lastYear>.png
#         output/species_trend_distributions/manifest_<firstYear>_<lastYear>.csv
#           (one row per species/model/scenario attempted)
#
# Resumable like 3d: a plot that already exists is skipped unless
# overwrite = TRUE. A species/model/scenario with fewer than
# min_routes_per_plot categorized routes is skipped (not enough data for a
# meaningful density) with a note in the manifest, not a fatal error.

library(here)
library(tidyverse)

here::i_am("6a_species_trend_distributions.R")

# Settings -------------------------------------------------------------------
firstYear <- 2010
lastYear  <- 2025

models    <- c("base", "anthro")
scenarios <- c(rcp45 = "RCP 4.5", rcp85 = "RCP 8.5")

# NULL = every species found in per_species_sdm/; otherwise a character
# vector of file-name species slugs (e.g. "Wood_Thrush"), handy for testing.
only_species <- NULL

overwrite <- FALSE

# Same rule as 3d/4c: drop routes whose own alpha_raw/beta_raw_space didn't converge.
require_route_converged <- TRUE

# Skip a species/model/scenario combo with fewer than this many categorized
# (Contraction/Stable/Expansion) routes -- not enough data for a meaningful
# density curve.
min_routes_per_plot <- 5

fig_width  <- 9
fig_height <- 6.5
fig_dpi    <- 200

sdm_dir <- here::here("output", "species_routes_covariates", "per_species_sdm")
out_dir <- here::here("output", "species_trend_distributions")

# Colours ---------------------------------------------------------------------
# Fill = exactly 3d_visualization_map.R's shift_cols, so this figure and the
# range-shift maps read as one system. Border/point/label/mean-triangle
# colors are manually darkened versions of the same hue -- Stable's map
# fill (#FFFF99) is a near-invisible pale yellow as a line/point/text color
# on a white background, so it needs a distinctly darker gold rather than a
# generic darkening formula.
shift_cols        <- c(Contraction = "#FF7F00", Stable = "#FFFF99", Expansion = "#19B2FF")
shift_border_cols <- c(Contraction = "#B35900", Stable = "#999900", Expansion = "#0F6E9C")
cat_levels <- c("Contraction", "Stable", "Expansion")

# Raster value -> range-shift class, same grouping as 3d/4c:
# Contraction = 1,2,3 | Stable = 4 | Expansion = 5,6,7 (0 = never suitable, excluded).
group_category <- function(x) {
  dplyr::case_when(
    x %in% c(1, 2, 3) ~ "Contraction",
    x == 4            ~ "Stable",
    x %in% c(5, 6, 7) ~ "Expansion",
    TRUE              ~ NA_character_
  )
}

# density() errors/degenerates with < 2 points or zero variance -- returns
# NULL in that case so callers can skip that curve rather than crash.
safe_density <- function(x) {
  if (length(x) < 2 || stats::var(x) == 0) return(NULL)
  tryCatch(stats::density(x), error = function(e) NULL)
}

# Load route trend tables ------------------------------------------------------
sdm_files <- list.files(sdm_dir, pattern = "_route_trends_sdm\\.csv$", full.names = TRUE)
if (length(sdm_files) == 0) {
  stop("No *_route_trends_sdm.csv files in ", sdm_dir, " -- run 3c_add_SDM_covariates.R first.")
}

trends <- map_dfr(sdm_files, function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  d$slug <- sub("_(base|anthro)_route_trends_sdm\\.csv$", "", basename(f))
  d
}) %>%
  filter(model %in% models, !is.na(trend), !is.na(trend_lci), !is.na(trend_uci))

if (require_route_converged) {
  n_before <- nrow(trends)
  trends <- trends %>% filter(route_converged %in% TRUE)
  cat("require_route_converged = TRUE -> dropped", n_before - nrow(trends), "of", n_before,
      "route/species/model rows.\n")
}

trends <- trends %>%
  mutate(
    trend_class = case_when(
      trend_lci > 0 ~ "Increase",
      trend_uci < 0 ~ "Decrease",
      TRUE          ~ "Not significant"
    ),
    credible = trend_class != "Not significant"
  )

species_slugs <- unique(trends$slug)
if (!is.null(only_species)) species_slugs <- intersect(species_slugs, only_species)
cat("Species to plot:", length(species_slugs), "| models:", paste(models, collapse = ", "),
    "| scenarios:", paste(scenarios, collapse = ", "), "\n\n")

# Plot builder -----------------------------------------------------------------
# df must already have one row per route with: trend, credible, category
# (Contraction/Stable/Expansion factor, NA rows already dropped).
plot_species_distribution <- function(df, sp_name, model, scenario_label) {
  dens_all <- safe_density(df$trend)

  cat_dfs <- lapply(cat_levels, function(lv) df %>% filter(category == lv))
  names(cat_dfs) <- cat_levels
  cat_dens <- lapply(cat_dfs, function(x) safe_density(x$trend))

  peak <- max(c(if (!is.null(dens_all)) dens_all$y else NA_real_,
               unlist(lapply(cat_dens, function(x) if (!is.null(x)) x$y else NA_real_))),
             na.rm = TRUE)
  x_max <- max(c(if (!is.null(dens_all)) dens_all$x else NA_real_,
                unlist(lapply(cat_dens, function(x) if (!is.null(x)) x$x else NA_real_)),
                df$trend),
              na.rm = TRUE)

  present <- cat_levels[sapply(cat_dfs, nrow) > 0]
  jitter_h <- peak * 0.075
  row_y <- setNames(-peak * 0.16 * seq_along(present), present)

  set.seed(1)  # reproducible jitter across re-runs
  for (lv in present) {
    n_lv <- nrow(cat_dfs[[lv]])
    cat_dfs[[lv]]$y_jit <- row_y[[lv]] + stats::runif(n_lv, -jitter_h, jitter_h)
  }

  mean_marks <- tibble(
    category = factor(present, levels = cat_levels),
    mean_trend = sapply(present, function(lv) mean(cat_dfs[[lv]]$trend))
  )

  n_counts <- sapply(cat_dfs, nrow)
  n_all <- nrow(df)

  p <- ggplot()

  # Category densities (Stable, then Contraction, then Expansion so the
  # narrower/more-saturated colors sit on top of the broader Stable fill).
  for (lv in intersect(c("Stable", "Contraction", "Expansion"), present)) {
    if (!is.null(cat_dens[[lv]])) {
      p <- p + geom_density(data = cat_dfs[[lv]], aes(x = trend), inherit.aes = FALSE,
                            fill = shift_cols[[lv]], color = shift_border_cols[[lv]],
                            alpha = if (lv == "Stable") 0.30 else 0.40, linewidth = 0.7)
    }
  }

  # Overall (every categorized route) drawn last, outline only, no fill --
  # a background reference silhouette that never competes for area.
  if (!is.null(dens_all)) {
    p <- p + geom_density(data = df, aes(x = trend), inherit.aes = FALSE,
                          fill = NA, color = "grey40", linewidth = 0.8, linetype = "22")
  }

  # Jittered per-route points + black credible ring, one row per category present.
  for (lv in present) {
    p <- p + geom_point(data = cat_dfs[[lv]], aes(x = trend, y = y_jit), inherit.aes = FALSE,
                        color = shift_cols[[lv]], size = 1.5, alpha = 0.75)
  }
  credible_pts <- bind_rows(cat_dfs[present]) %>% filter(credible)
  if (nrow(credible_pts) > 0) {
    p <- p + geom_point(data = credible_pts, aes(x = trend, y = y_jit), inherit.aes = FALSE,
                        shape = 1, size = 2.5, stroke = 0.6, color = "black")
  }

  # Mean-trend triangle per category, at the density baseline (y = 0).
  p <- p +
    geom_point(data = mean_marks, aes(x = mean_trend, y = 0, fill = category), inherit.aes = FALSE,
               shape = 24, size = 3.6, color = "black", stroke = 0.5) +
    scale_fill_manual(values = shift_cols, guide = "none") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey20", linewidth = 0.6)

  # Right-aligned "<Category> (n=..)" label at each category's own row.
  for (lv in present) {
    p <- p + annotate("text", x = x_max * 0.98, y = row_y[[lv]],
                      label = paste0(lv, " (n=", n_counts[[lv]], ")"), hjust = 1, size = 3.6,
                      color = shift_border_cols[[lv]], fontface = "bold")
  }

  count_str <- paste0(n_all, " routes total (",
                      paste0(n_counts[present], " ", tolower(present), collapse = ", "), ")")

  p <- p +
    labs(
      title = sp_name,
      subtitle = paste0(scenario_label, " range-shift categories | ", model, " model | ", count_str),
      x = "Route trend (annual % change)",
      y = "Density",
      caption = paste0("Black ring = route's own 90% credible interval excludes zero  |  ",
                       "▲ = category mean trend  |  dashed grey outline = all categorized routes (background)")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold", size = 16),
      plot.subtitle    = element_text(size = 10.5, color = "grey30"),
      plot.caption     = element_text(size = 8.5, color = "grey40"),
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      panel.grid.minor = element_blank()
    )

  p
}

# Main loop --------------------------------------------------------------------
manifest <- list()
add_manifest <- function(slug, model, scenario, n_routes, n_contraction, n_stable, n_expansion, status) {
  manifest[[length(manifest) + 1]] <<- data.frame(
    species = slug, model = model, scenario = scenario, n_routes = n_routes,
    n_contraction = n_contraction, n_stable = n_stable, n_expansion = n_expansion, status = status
  )
}

for (sp in species_slugs) {
  sp_all <- trends %>% filter(slug == sp)
  sp_name <- unique(sp_all$species)[1]

  jobs <- expand.grid(model = intersect(models, unique(sp_all$model)),
                      scenario = names(scenarios), stringsAsFactors = FALSE) %>%
    mutate(file = file.path(out_dir, model, scenario,
                            paste0(sp, "_", model, "_", scenario, "_", firstYear, "_", lastYear, ".png")))
  todo <- if (overwrite) jobs else jobs %>% filter(!file.exists(file))
  if (nrow(todo) == 0) { cat("Skipping (all plots exist):", sp, "\n"); next }

  cat("Processing:", sp_name, "\n")

  for (i in seq_len(nrow(todo))) {
    m  <- todo$model[i]
    sc <- todo$scenario[i]
    out_file <- todo$file[i]

    df <- sp_all %>%
      filter(model == m, !is.na(.data[[sc]])) %>%
      mutate(category = factor(group_category(.data[[sc]]), levels = cat_levels)) %>%
      filter(!is.na(category))

    n_counts <- table(factor(df$category, levels = cat_levels))

    if (nrow(df) < min_routes_per_plot) {
      add_manifest(sp, m, sc, nrow(df), n_counts[["Contraction"]], n_counts[["Stable"]],
                  n_counts[["Expansion"]], "too_few_routes")
      next
    }

    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

    status <- tryCatch({
      p <- plot_species_distribution(df, sp_name, m, scenarios[[sc]])
      ggsave(out_file, p, width = fig_width, height = fig_height, dpi = fig_dpi, device = ragg::agg_png)
      "ok"
    }, error = function(e) {
      message("WARNING: ", sp, " ", m, " ", sc, " -- plot failed: ", conditionMessage(e))
      "error"
    })
    add_manifest(sp, m, sc, nrow(df), n_counts[["Contraction"]], n_counts[["Stable"]],
                n_counts[["Expansion"]], status)
    if (status == "ok") cat("  Wrote:", basename(out_file), "\n")
  }
}

# Summary ------------------------------------------------------------------
if (length(manifest) > 0) {
  manifest_df <- bind_rows(manifest)
  manifest_csv <- file.path(out_dir, paste0("manifest_", firstYear, "_", lastYear, ".csv"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(manifest_df, manifest_csv, row.names = FALSE)
  cat("\n=== Species trend distributions done ===\n")
  print(table(manifest_df$status))
  cat("Manifest:", manifest_csv, "\n")
} else {
  cat("\nNothing to do -- every requested plot already exists (set overwrite <- TRUE to redraw).\n")
}
