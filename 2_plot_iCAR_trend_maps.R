## =============================================================================
## Plot iCAR route-level trend maps
## Time period: 2010-2024
##
## Reads the per-species route-trend CSVs produced by
## 1_species_iCAR_2010_2024.R and generates one trend map per species.
##
## Input  : output/species_routes/<group>_<species>_route_trends.csv
##          columns: species, species_code, group, route, routeF,
##                   latitude, longitude, trend, trend_lci, trend_uci,
##                   rel_abundance
## Output : figures/<species>_iCAR_trend_map_<firstYear>_<lastYear>.png (+ .pdf)
## =============================================================================

library(bbsBayes2)
library(tidyverse)
library(sf)
library(here)

here::i_am("2_plot_iCAR_trend_maps.R")

# Settings (match 1_species_iCAR_2010_2024.R) -----------------------------
land_cover <- "grasslands"
firstYear  <- 2010
lastYear   <- 2024
strat      <- "bbs_usgs"

# Also write PDF alongside PNG?
write_pdf <- TRUE

# Directories
sp_out_dir <- here::here("output", "species_routes")
fig_dir    <- here::here("figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Helper: convert species name to file-safe format ------------------------
species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

# Trend categories / palette ----------------------------------------------
trend_breaks  <- c(-Inf, -4, -2, -0.5, 0.5, 2, 4, Inf)
trend_labels  <- c("< -4%", "-4 to -2%", "-2 to -0.5%", "-0.5 to 0.5%",
                   "0.5 to 2%", "2 to 4%", "> 4%")
trend_colours <- c("#d73027", "#fc8d59", "#fee08b", "#ffffbf",
                   "#d9ef8b", "#91cf60", "#1a9850")

# Base strata map (shared across species) ---------------------------------
strata_map <- load_map(strat)

# Locate per-species CSVs --------------------------------------------------
csv_files <- list.files(sp_out_dir,
                        pattern = paste0("^", land_cover, "_.*_route_trends\\.csv$"),
                        full.names = TRUE)

if (length(csv_files) == 0) {
  stop("No route-trend CSVs found in ", sp_out_dir,
       " — run 1_species_iCAR_2010_2024.R first.")
}

cat("=== Plotting iCAR trend maps ===\n")
cat("Group:", land_cover, " | Period:", firstYear, "-", lastYear, "\n")
cat("CSVs found:", length(csv_files), "\n")

# ==========================================================================
# Loop: one trend map per species CSV
# ==========================================================================
for (csv in csv_files) {

  route_trends <- read.csv(csv, stringsAsFactors = FALSE)
  if (nrow(route_trends) == 0) {
    cat("  Skipping (empty):", basename(csv), "\n")
    next
  }

  sp   <- route_trends$species[1]
  sp_f <- species_to_f(sp)

  cat("\n  ", sp, "(", nrow(route_trends), "routes )\n")

  # Build sf points from lat/lon, projected to the strata-map CRS
  trend_sf <- st_as_sf(route_trends, coords = c("longitude", "latitude"),
                       crs = 4326) %>%
    st_transform(crs = st_crs(strata_map))

  trend_sf$trend_cat <- cut(trend_sf$trend,
                            breaks = trend_breaks, labels = trend_labels)

  # Background: only strata that contain routes for this species
  realized_strata_map <- strata_map %>%
    filter(strata_name %in% unique(route_trends$strata_name))
  if (nrow(realized_strata_map) == 0) {
    # CSV may not carry strata_name; fall back to clipping near the points
    realized_strata_map <- strata_map[st_buffer(st_union(trend_sf), 50000), ]
  }

  p_trend <- ggplot() +
    geom_sf(data = realized_strata_map, fill = grey(0.95), colour = grey(0.8)) +
    geom_sf(data = trend_sf,
            aes(colour = trend_cat, size = rel_abundance), alpha = 0.8) +
    scale_colour_manual(values = trend_colours, name = "Annual trend", drop = FALSE) +
    scale_size_continuous(range = c(0.5, 4), name = "Relative\nabundance") +
    labs(title = paste(sp, "— iCAR route-level trends", firstYear, "-", lastYear),
         subtitle = paste("Routes:", nrow(route_trends),
                          " | Median trend:", round(median(route_trends$trend), 2), "%")) +
    theme_bw() +
    theme(legend.position = "right",
          text = element_text(family = "serif", size = 11))

  png(file.path(fig_dir, paste0(sp_f, "_iCAR_trend_map_",
                                firstYear, "_", lastYear, ".png")),
      width = 11, height = 8.5, units = "in", res = 300)
  print(p_trend)
  dev.off()

  if (write_pdf) {
    pdf(file.path(fig_dir, paste0(sp_f, "_iCAR_trend_map_",
                                  firstYear, "_", lastYear, ".pdf")),
        width = 11, height = 8.5)
    print(p_trend)
    dev.off()
  }

  cat("    Map saved to figures/\n")
}

cat("\n=== Done ===\n")
cat("Trend maps written to:", fig_dir, "\n")
