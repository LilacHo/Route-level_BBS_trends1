# 6_visualization_publishable.R
#
library(here)
library(tidyverse)
library(sf)

here::i_am("6_visualization_publishable.R")

# Settings -------------------------------------------------------------------
firstYear <- 2010
lastYear  <- 2025
model_tag <- "anthro"   # which fit's route set to map ("base" or "anthro")

# Same rule as 3d/4c: only count route x species rows whose route-level
# parameters converged.
require_route_converged <- TRUE

# Drop US-RIIS-listed (non-native) species.
drop_invasive <- TRUE

# Colour scale for the average route trend is symmetric around 0 and capped at
# this quantile of |mean trend|; routes beyond it get the end colour.
trend_cap_quantile <- 0.98

fig_width <- 10
fig_dpi   <- 300

in_csv  <- here::here("output", "species_routes_covariates",
                      paste0("all_route_trends_all_species_", model_tag, "_", firstYear, "_", lastYear, ".csv"))
out_dir <- here::here("output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Base map (identical to 3d) -------------------------------------------------
map_crs <- "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

states_map <- bbsBayes2::load_map("prov_state") %>% st_transform(map_crs)
us_outline <- states_map %>%
  filter(country_code == "US") %>%
  st_make_valid() %>%
  st_union() %>%
  st_make_valid()
# North America only: other continents' polygons wrap badly in this projection
# and covered the ocean.
world_map <- spData::world %>%
  filter(continent == "North America", !iso_a2 %in% c("US", "CA")) %>%
  st_transform(map_crs)

# Routes ---------------------------------------------------------------------
trends <- read_csv(in_csv, show_col_types = FALSE) %>%
  filter(model == model_tag, !is.na(route), !is.na(latitude), !is.na(longitude))
if (require_route_converged) trends <- trends %>% filter(route_converged %in% TRUE)

if (drop_invasive) {
  traits <- read_csv(here::here("data", "spp_names_codes_group_aou_taxa_redlist_diet_migration_invasive.csv"),
                     show_col_types = FALSE)
  non_native <- traits %>% filter(!is.na(invasive), invasive != "Not listed") %>% pull(Code)
  dropped <- trends %>% filter(species_code %in% non_native) %>% distinct(species) %>% pull(species)
  cat("Dropping", length(dropped), "non-native species:", paste(sort(dropped), collapse = ", "), "\n")
  trends <- trends %>% filter(!species_code %in% non_native)
}

routes <- trends %>%
  group_by(route) %>%
  summarise(latitude  = first(latitude),
            longitude = first(longitude),
            n_species = n_distinct(species_code),
            mean_trend = mean(trend),
            .groups = "drop")

routes_sf <- st_as_sf(routes, coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(map_crs)

cat(nrow(routes), "routes;", n_distinct(trends$species_code), "species;",
    "species per route: median", median(routes$n_species),
    "range", min(routes$n_species), "-", max(routes$n_species), "\n")

# Extent: the routes plus a small margin
bb  <- st_bbox(routes_sf)
pad <- 0.03 * max(bb[["xmax"]] - bb[["xmin"]], bb[["ymax"]] - bb[["ymin"]])
xlim <- c(bb[["xmin"]] - pad, bb[["xmax"]] + pad)
ylim <- c(bb[["ymin"]] - pad, bb[["ymax"]] + pad)
fig_height <- fig_width * diff(ylim) / diff(xlim) + 0.6  # + room for title (legend sits inside the map)

base_layers <- list(
  geom_sf(data = world_map,  fill = "grey95", colour = "grey75", linewidth = 0.2),
  geom_sf(data = states_map, fill = "grey95", colour = "grey60", linewidth = 0.2),
  geom_sf(data = us_outline, fill = NA,       colour = "grey35", linewidth = 0.4),
  theme_bw(base_size = 12),
  theme(panel.background = element_rect(fill = "aliceblue"),
        panel.grid = element_line(colour = "white", linewidth = 0.2),
        axis.title = element_blank())
)

# Added last in each plot: a geom_sf() added after a coord_sf() replaces it
# with a default (unzoomed) one. default_crs = map_crs so xlim/ylim are read
# as Albers metres, not lon/lat.
map_coord <- coord_sf(crs = map_crs, default_crs = map_crs, xlim = xlim, ylim = ylim, expand = FALSE)

save_fig <- function(p, stem) {
  p <- p + map_coord
  ggsave(file.path(out_dir, paste0(stem, ".png")),
         p, width = fig_width, height = fig_height, dpi = fig_dpi)
}

# Fig 1: routes, size = species, colour = average trend ---------------------
trend_cap <- ceiling(quantile(abs(routes$mean_trend), trend_cap_quantile))
fig1 <- ggplot() +
  base_layers +
  # big points first so small ones stay visible on top
  geom_sf(data = arrange(routes_sf, desc(n_species)),
          aes(size = n_species, fill = mean_trend),
          shape = 21, colour = "grey25", stroke = 0.15, alpha = 0.9) +
  scale_fill_distiller(palette = "RdBu", direction = 1,
                       limits = c(-trend_cap, trend_cap), oob = scales::squish,
                       name = "Average route trend (% per year)",
                       guide = guide_colourbar(barwidth = 10, barheight = 0.6, title.position = "top",
                                               order = 1)) +
  scale_size(range = c(0.4, 3.2), name = "Number of species modelled per route",
             guide = "none") +
             # guide = guide_legend(title.position = "top", order = 2,
             #                      override.aes = list(fill = "grey70"))) +
  labs(
    # title = "North American Breeding Bird Survey routes",
       subtitle = paste0(format(nrow(routes), big.mark = ","), " routes, ",
                         n_distinct(trends$species_code), " native species, ", firstYear, "-", lastYear#,
                         # if (drop_invasive) " (non-native species excluded)" else ""
                         )) +
  # Legend inside the map, bottom-left (Pacific / Mexico: no routes there)
  theme(legend.position = "inside",
        legend.position.inside = c(0.02, 0.02),
        legend.justification = c(0, 0),
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.box.just = "left",
        legend.spacing.y = unit(0.2, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.box.background = element_rect(fill = alpha("white", 0.85), colour = "grey60", linewidth = 0.3),
        legend.box.margin = margin(4, 6, 4, 6),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))
save_fig(fig1, "figure1")
