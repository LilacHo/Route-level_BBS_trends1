# 3d_visualization_map.R
#
# One map per species x model (base / anthro) x climate scenario (RCP 4.5 /
# RCP 8.5) showing, on a single panel:
#   - the SDM range shift (filled area): Contraction / Stable / Expansion,
#     from the 2025 classified-change raster (Bateman et al. 2020), and
#   - each BBS route's population trend (points): colour = Decrease /
#     Increase / Not significant, size = |annual % change|.
#
# Replaces the earlier single-species (Blue Jay) MARSS version of this script;
# the map layout is the same idea, but the trends now come from this
# project's iCAR route-level fits (2c) and the SDM category from the rasters
# 3c already extracts against, so every species/model/scenario is drawn the
# same way with no per-species code.
#
# Range-shift classes (raster values -> map category):
#   1-3 (extirpation, worsening, slightly worsening) -> Contraction
#   4   (neutral)                                    -> Stable
#   5-7 (slightly improving, improving, colonization)-> Expansion
#   0   (never suitable)                             -> not drawn
# This is the same Contraction / Stable / Expansion grouping 4c uses in its
# PART 3-4.
#
# Trend class (per route, from the 90% credible interval 2c writes as
# trend_lci / trend_uci):
#   trend_lci > 0 -> Increase | trend_uci < 0 -> Decrease | otherwise Not significant
#
# Reads:  output/species_routes_covariates/per_species_sdm/*_route_trends_sdm.csv
#           (3c; species, group, route coords, trend, trend_lci/uci, route_converged)
#         data/rcp{45,85}_<group>/<code>/<group>_<code>_breeding_2025_{45,85}_ENSEMBLE_classifiedchange.tif
#         bbsBayes2's "prov_state" map + spData::world for the base map
# Writes: output/species_maps/<model>/<scenario>/<species>_<model>_<scenario>_<firstYear>_<lastYear>.png
#         output/species_maps/map_manifest_<firstYear>_<lastYear>.csv (one row per map attempted)
#
# Resumable like 1c: a map that already exists is skipped unless
# overwrite = TRUE, so this can be stopped and restarted. A species whose
# raster is missing for a scenario is skipped for that scenario with a
# warning (same convention as 3c), not treated as fatal.

library(here)
library(tidyverse)
library(terra)
library(sf)

here::i_am("3d_visualization_map.R")

# Settings -------------------------------------------------------------------
firstYear <- 2010
lastYear  <- 2025

models    <- c("base", "anthro")
scenarios <- c(rcp45 = "RCP 4.5", rcp85 = "RCP 8.5")  # names = raster folder prefix

# NULL = every species found in per_species_sdm/; otherwise a character vector
# of file-name species slugs (e.g. "Eastern_Kingbird"), handy for testing.
only_species <- NULL

overwrite <- FALSE

# Same rule as 4c: drop routes whose own alpha_raw/beta_raw_space didn't converge.
require_route_converged <- TRUE

# Point size is |annual % change| squished into [0, size_limit] so sizes are
# comparable between maps (a species with huge trends can't shrink everyone
# else's points to nothing).
size_limit <- 10

# Range raster is 1 km; it's downsampled (modal class) to about this many
# columns across the plotted extent -- plenty for a 10-inch, 300-dpi figure
# and keeps the dissolved range polygons light.
target_cols <- 2200

# Zoom: the map never shows more than the contiguous (lower-48) US, and zooms
# in further to just the routes + range that fall inside it (so a Southwest
# species gets a Southwest map, not a whole-CONUS map). Canada/Alaska/Mexico
# parts of a range and routes outside the lower 48 are cropped out of the
# view. A species with no routes and no range inside the lower 48 (e.g. an
# Alaska-only bird) falls back to the whole-range extent below.
conus_zoom <- TRUE

# Fallback extent (conus_zoom = FALSE, or nothing inside the lower 48): the
# routes plus the central (1 - 2 * extent_trim) of the range cells'
# coordinates, so scattered outlier patches far from the main range don't
# stretch the map into empty space. 0 = use the full range.
extent_trim <- 0.01

fig_width <- 10   # inches; height follows the plotted extent's aspect ratio
fig_dpi   <- 300

sdm_dir <- here::here("output", "species_routes_covariates", "per_species_sdm")
out_dir <- here::here("output", "species_maps")

# Colours (range shift as in the Blue Jay reference map) --------------------
shift_cols <- c("Contraction" = "#FF7F00", "Stable" = "#FFFF99", "Expansion" = "#19B2FF")
trend_cols <- c("Decrease" = "#E51932", "Increase" = "#654CFF", "Not significant" = "#8FD3C0")

# Raster value -> range-shift class (1 = Contraction, 2 = Stable, 3 = Expansion);
# unlisted values (0 = never suitable) become NA and are not drawn.
shift_rcl <- cbind(is = 1:7, becomes = c(1, 1, 1, 2, 3, 3, 3))
shift_labels <- names(shift_cols)

# Base map -------------------------------------------------------------------
# The SDM rasters are all in this Albers equal-area conic, so the whole map
# is drawn in it and nothing needs reprojecting.
map_crs <- "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

states_map <- bbsBayes2::load_map("prov_state") %>% st_transform(map_crs)
# Dissolved US outline, drawn thicker than the state/province lines.
us_outline <- states_map %>%
  filter(country_code == "US") %>%
  st_make_valid() %>%
  st_union() %>%
  st_make_valid()
world_map  <- spData::world %>%
  st_transform(map_crs) %>%
  filter(!iso_a2 %in% c("US", "CA"))   # US/Canada are drawn from states_map

# Lower-48 bounding box (map CRS) that conus_zoom clips the view to.
conus_bbox <- spData::us_states %>% st_transform(map_crs) %>% st_bbox()
conus_ext  <- terra::ext(conus_bbox[["xmin"]], conus_bbox[["xmax"]], conus_bbox[["ymin"]], conus_bbox[["ymax"]])

# Load route trend tables -----------------------------------------------------
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

# A handful of rows across the pipeline carry no route/lat/lon (an upstream
# join gap -- e.g. a route dropped from route_keys but not from the trend
# table); such a row can't be placed on a map and crashes st_as_sf() below,
# so it's dropped here with a warning rather than failing the whole run.
missing_coords <- is.na(trends$latitude) | is.na(trends$longitude)
if (any(missing_coords)) {
  message("WARNING: dropped ", sum(missing_coords), " row(s) with missing latitude/longitude ",
          "(species: ", paste(sort(unique(trends$slug[missing_coords])), collapse = ", "), ").")
  trends <- trends %>% filter(!missing_coords)
}

if (require_route_converged) {
  n_before <- nrow(trends)
  trends <- trends %>% filter(route_converged %in% TRUE)
  cat("require_route_converged = TRUE -> dropped", n_before - nrow(trends), "of", n_before,
      "route/species/model rows.\n")
}

trends <- trends %>%
  mutate(trend_class = case_when(
    trend_lci > 0 ~ "Increase",
    trend_uci < 0 ~ "Decrease",
    TRUE          ~ "Not significant"
  ),
  trend_class = factor(trend_class, levels = names(trend_cols)))

species_slugs <- unique(trends$slug)
if (!is.null(only_species)) species_slugs <- intersect(species_slugs, only_species)
cat("Species to map:", length(species_slugs), "| models:", paste(models, collapse = ", "),
    "| scenarios:", paste(scenarios, collapse = ", "), "\n\n")

# Helpers --------------------------------------------------------------------

# Range-shift raster for one species/scenario as dissolved polygons (sf),
# cropped to the plotted extent. `route_bbox` = bbox of the routes to frame
# (those inside the lower 48 when conus_zoom, else all of them), or NULL if
# there are none; `all_route_bbox` = bbox of every route, used by the
# fallback extent. Returns list(shapes, xlim, ylim), or NULL if the raster is
# missing.
prepare_range <- function(tif, route_bbox, all_route_bbox) {
  if (!file.exists(tif)) return(NULL)

  r <- rast(tif)

  # A handful of rasters (confirmed: some rcp45_waterbirds and
  # rcp45_western_forests files) encode the 0-7 classes as 0/10000/.../70000
  # instead -- an upstream data bug, not a different legend. shift_rcl only
  # matches 1-7, so left as-is these rasters classify to all-NA and silently
  # draw no range at all. Detected via the actual max value (minmax()'s
  # cached header stats can be missing -- global() forces a real scan when
  # that happens) and corrected by rescaling back to the 0-7 legend.
  rmax <- terra::minmax(r)[2, 1]
  if (is.na(rmax)) rmax <- terra::global(r, "max", na.rm = TRUE)[1, 1]
  if (!is.na(rmax) && rmax > 7) {
    message("NOTE: ", basename(tif), " -- raster values scaled x10000 (max = ", rmax,
            "); rescaling back to the 0-7 legend.")
    r <- round(r / 10000)
  }

  if (!terra::same.crs(r, map_crs)) r <- terra::project(r, map_crs, method = "near")

  # Coarse (10 km) range cells, just to find where the range is.
  range_cells <- function(x) {
    cls <- terra::classify(x, shift_rcl, others = NA)
    coarse <- terra::aggregate(cls, fact = 10, fun = "modal", na.rm = TRUE)
    terra::as.data.frame(coarse, xy = TRUE, na.rm = TRUE)
  }

  use_conus <- conus_zoom
  if (use_conus) {
    cc <- range_cells(terra::crop(r, conus_ext))
    if (is.null(route_bbox) && nrow(cc) == 0) use_conus <- FALSE
  }
  if (!use_conus) {
    cc <- range_cells(r)
    route_bbox <- all_route_bbox
  }

  # Fallback extent uses quantiles of the range cells' coordinates
  # (extent_trim), not the strict min/max, so isolated patches far from the
  # main range (e.g. stray cells in Alaska) can't stretch the map into empty
  # space. Inside the lower 48 the range is already bounded, so the strict
  # min/max is used there.
  trim <- if (use_conus) 0 else extent_trim

  xmin <- route_bbox[["xmin"]]; xmax <- route_bbox[["xmax"]]
  ymin <- route_bbox[["ymin"]]; ymax <- route_bbox[["ymax"]]
  if (nrow(cc) > 0) {
    qx <- quantile(cc$x, c(trim, 1 - trim), names = FALSE)
    qy <- quantile(cc$y, c(trim, 1 - trim), names = FALSE)
    xmin <- min(c(xmin, qx[1])); xmax <- max(c(xmax, qx[2]))
    ymin <- min(c(ymin, qy[1])); ymax <- max(c(ymax, qy[2]))
  }
  padx <- (xmax - xmin) * 0.04; pady <- (ymax - ymin) * 0.04
  plot_ext <- terra::ext(xmin - padx, xmax + padx, ymin - pady, ymax + pady)

  # Crop the drawn range from the FULL raster (not the lower-48-cropped one
  # used only to find the extent), so the range runs to the panel edge
  # instead of ending in a straight cut at the lower-48 bounding box.
  cropped <- terra::classify(terra::crop(r, plot_ext), shift_rcl, others = NA)
  fact <- max(1, floor(ncol(cropped) / target_cols))
  if (fact > 1) cropped <- terra::aggregate(cropped, fact = fact, fun = "modal", na.rm = TRUE)

  # Dissolve the classes into polygons (one multipolygon per class) rather
  # than drawing cells: crisp edges, and geom_raster under coord_sf warns
  # about "uneven intervals" for these sparse grids.
  names(cropped) <- "shift"
  shapes <- sf::st_as_sf(terra::as.polygons(cropped, dissolve = TRUE, na.rm = TRUE))
  shapes$shift <- factor(shapes$shift, levels = 1:3, labels = shift_labels)

  list(shapes = shapes, xlim = c(plot_ext$xmin, plot_ext$xmax), ylim = c(plot_ext$ymin, plot_ext$ymax))
}

plot_species_map <- function(range, routes_sf, species, model, scenario_label) {
  # Only routes inside the plotted window are counted in the caption.
  xy <- sf::st_coordinates(routes_sf)
  in_view <- xy[, 1] >= range$xlim[1] & xy[, 1] <= range$xlim[2] &
             xy[, 2] >= range$ylim[1] & xy[, 2] <= range$ylim[2]
  n_class <- table(routes_sf$trend_class[in_view])
  n_outside <- sum(!in_view)

  # Invisible points, one per trend class, so the legend always lists all
  # three classes even when a species has no routes in one of them (the
  # legend override below makes their keys visible).
  legend_dummy <- data.frame(x = mean(range$xlim), y = mean(range$ylim),
                             trend_class = factor(names(trend_cols), levels = names(trend_cols)))

  p <- ggplot() +
    geom_sf(data = world_map,  fill = NA, color = "grey90", linewidth = 0.3) +
    geom_sf(data = range$shapes, aes(fill = shift), color = NA) +
    geom_sf(data = states_map, fill = NA, color = "#BFBFBF", linewidth = 0.25) +
    geom_sf(data = us_outline, fill = NA, color = "grey35", linewidth = 0.8) +
    geom_point(data = legend_dummy, aes(x = x, y = y, color = trend_class), alpha = 0) +
    geom_sf(data = routes_sf, aes(color = trend_class, size = pmin(abs(trend), size_limit)), alpha = 0.7) +
    scale_fill_manual(name = "Range Shift", values = shift_cols, drop = FALSE) +
    scale_color_manual(name = "Population Trend", values = trend_cols, drop = FALSE) +
    scale_size_continuous(range = c(0.8, 4), limits = c(0, size_limit), guide = "none") +
    guides(fill  = guide_legend(order = 1),
           color = guide_legend(order = 2, override.aes = list(size = 3, alpha = 1))) +
    coord_sf(xlim = range$xlim, ylim = range$ylim, expand = FALSE) +
    labs(title = species,
         subtitle = paste0(scenario_label, " range shift (2025); BBS route trends ", firstYear, "-", lastYear,
                           " | model: ", model),
         caption = paste0(sum(in_view), " routes (", n_class[["Decrease"]], " decreasing, ",
                          n_class[["Increase"]], " increasing, ", n_class[["Not significant"]],
                          " not significant; 90% CI)",
                          if (n_outside > 0) paste0("; ", n_outside, " more outside the map") else "")) +
    theme_light(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold", size = 17),
      plot.subtitle    = element_text(size = 11, color = "grey30"),
      plot.caption     = element_text(size = 9, color = "grey40"),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      axis.title       = element_blank(),
      legend.position        = "inside",
      legend.position.inside = c(0.02, 0.02),
      legend.justification   = c(0, 0),
      legend.box             = "vertical",
      legend.box.just        = "left",
      legend.background      = element_blank(),
      legend.key             = element_blank(),
      legend.title           = element_text(face = "bold", size = 12),
      legend.text            = element_text(size = 11)
    )

  aspect <- diff(range$ylim) / diff(range$xlim)
  list(plot = p, height = min(max(fig_width * aspect + 1.4, 5), 11))
}

# Main loop ------------------------------------------------------------------
manifest <- list()
add_manifest <- function(slug, model, scenario, n_routes, status) {
  manifest[[length(manifest) + 1]] <<- data.frame(species = slug, model = model, scenario = scenario,
                                                  n_routes = n_routes, status = status)
}

for (sp in species_slugs) {
  sp_all <- trends %>% filter(slug == sp)
  code    <- unique(sp_all$species_code)[1]
  group   <- unique(sp_all$group)[1]
  sp_name <- unique(sp_all$species)[1]

  jobs <- expand.grid(model = intersect(models, unique(sp_all$model)),
                      scenario = names(scenarios), stringsAsFactors = FALSE) %>%
    mutate(file = file.path(out_dir, model, scenario,
                            paste0(sp, "_", model, "_", scenario, "_", firstYear, "_", lastYear, ".png")))
  todo <- if (overwrite) jobs else jobs %>% filter(!file.exists(file))
  if (nrow(todo) == 0) { cat("Skipping (all maps exist):", sp, "\n"); next }

  cat("Processing:", sp_name, "[", group, "]\n")

  all_routes_sf <- st_as_sf(sp_all, coords = c("longitude", "latitude"), crs = 4326) %>%
    st_transform(map_crs)
  all_route_bbox <- st_bbox(all_routes_sf)
  xy_all <- st_coordinates(all_routes_sf)
  in_conus <- xy_all[, 1] >= conus_bbox[["xmin"]] & xy_all[, 1] <= conus_bbox[["xmax"]] &
              xy_all[, 2] >= conus_bbox[["ymin"]] & xy_all[, 2] <= conus_bbox[["ymax"]]
  route_bbox <- if (conus_zoom) {
    if (any(in_conus)) st_bbox(all_routes_sf[in_conus, ]) else NULL
  } else all_route_bbox

  for (sc in unique(todo$scenario)) {
    suffix <- sub("^rcp", "", sc)
    tif <- here::here("data", paste0(sc, "_", group), code,
                      paste0(group, "_", code, "_breeding_2025_", suffix, "_ENSEMBLE_classifiedchange.tif"))

    range <- tryCatch(prepare_range(tif, route_bbox, all_route_bbox), error = function(e) {
      message("WARNING: ", sp, " ", sc, " -- raster prep failed: ", conditionMessage(e)); NULL
    })
    if (is.null(range)) {
      message("WARNING: ", sp, " ", sc, " -- raster missing/unreadable (", tif, ") -- skipping.")
      for (m in todo$model[todo$scenario == sc]) add_manifest(sp, m, sc, NA, "no_raster")
      next
    }

    for (m in todo$model[todo$scenario == sc]) {
      routes_sf <- all_routes_sf %>%
        filter(model == m) %>%
        arrange(trend_class != "Not significant")  # significant points drawn on top
      out_file <- todo$file[todo$model == m & todo$scenario == sc]
      dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

      status <- tryCatch({
        pm <- plot_species_map(range, routes_sf, sp_name, m, scenarios[[sc]])
        ggsave(out_file, pm$plot, width = fig_width, height = pm$height, dpi = fig_dpi,
               device = ragg::agg_png)
        "ok"
      }, error = function(e) {
        message("WARNING: ", sp, " ", m, " ", sc, " -- plot failed: ", conditionMessage(e))
        "error"
      })
      add_manifest(sp, m, sc, nrow(routes_sf), status)
      if (status == "ok") cat("  Wrote:", basename(out_file), "\n")
    }
  }
}

# Summary --------------------------------------------------------------------
if (length(manifest) > 0) {
  manifest_df <- bind_rows(manifest)
  manifest_csv <- file.path(out_dir, paste0("map_manifest_", firstYear, "_", lastYear, ".csv"))
  write.csv(manifest_df, manifest_csv, row.names = FALSE)
  cat("\n=== Species maps done ===\n")
  print(table(manifest_df$status))
  cat("Manifest:", manifest_csv, "\n")
} else {
  cat("\nNothing to do -- every requested map already exists (set overwrite <- TRUE to redraw).\n")
}
