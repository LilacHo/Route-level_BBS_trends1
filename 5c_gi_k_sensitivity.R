# 5c_gi_k_sensitivity.R
#
# Sensitivity of 5b_refugia_hotspot_stats.R's Getis-Ord Gi* hot spots to the
# neighborhood size k (number of nearest-neighbor routes). k sets the spatial
# scale of "local", and 5b's k = 8 was a convention, not a tuned choice. This
# reruns the same Gi* (same inputs, same weights construction as 5b) at
# several k and asks which hot/cold spots persist across all of them --
# clusters that survive every k are robust to the scale choice; ones that
# appear at only one k are scale artifacts.
#
# Reads:  output/species_routes_covariates/hotspots/refugia_stats/*_route_gi.csv
#         (written by 5b_refugia_hotspot_stats.R -- run that first)
# Writes: output/species_routes_covariates/hotspots/refugia_stats/k_sensitivity/
#           gi_k_sensitivity_summary.csv          (one row per combo x attribute x k)
#           gi_k_sensitivity_route_persistence.csv (one row per route x combo x attribute)
#           plots/<combo>_k_persistence_map.png    (rate | magnitude panels)

library(here)
library(tidyverse)
library(sf)
library(spdep)

here::i_am("5c_gi_k_sensitivity.R")

# Settings -------------------------------------------------------------------
k_values    <- c(5, 8, 12, 20, 30)
k_reference <- 8      # 5b's k
z_hot       <- 1.96   # |z| >= 1.96 (95%) counts as hot / cold for persistence
z_hot99     <- 2.58
fdr_q       <- 0.05   # Benjamini-Hochberg threshold, same as 5b
min_n_category <- 3   # same as 5b (rate routes need >= this many predicted species)

in_dir  <- here::here("output", "species_routes_covariates", "hotspots", "refugia_stats")
out_dir <- file.path(in_dir, "k_sensitivity")
plot_dir <- file.path(out_dir, "plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data("us_states", package = "spData")
us_states <- st_transform(us_states, 4326)
us_bbox <- st_bbox(us_states)
map_xlim <- c(us_bbox["xmin"] - 1, us_bbox["xmax"] + 1)
map_ylim <- c(us_bbox["ymin"] - 1, us_bbox["ymax"] + 1)

hav_km <- function(lon1, lat1, lon2, lat2) {
  r <- 6371; p <- pi / 180
  a <- sin((lat2 - lat1) * p / 2)^2 +
    cos(lat1 * p) * cos(lat2 * p) * sin((lon2 - lon1) * p / 2)^2
  2 * r * asin(pmin(1, sqrt(a)))
}

# Same construction as 5b's run_getis_ord(): planar k-nearest neighbors on
# lon/lat, self included (Gi*), binary weights. Also returns each route's
# distance (km) to its k-th nearest neighbor -- the actual spatial scale.
gi_star_at_k <- function(df, value_col, k) {
  coords <- as.matrix(df[, c("longitude", "latitude")])
  nb <- suppressWarnings(knn2nb(knearneigh(coords, k = k)))
  kth_km <- vapply(seq_along(nb), function(i) {
    max(hav_km(coords[i, 1], coords[i, 2], coords[nb[[i]], 1], coords[nb[[i]], 2]))
  }, numeric(1))
  z <- as.numeric(localG(df[[value_col]], nb2listw(include.self(nb), style = "B")))
  # Benjamini-Hochberg FDR across all routes tested at this k, as in 5b.
  q <- p.adjust(2 * pnorm(-abs(z)), method = "BH")
  list(z = z, q = q, kth_km = kth_km)
}

files <- list.files(in_dir, pattern = "_route_gi\\.csv$", full.names = TRUE)
# Expansion section is paused in 5b; skip any leftover expansion outputs from
# earlier runs. Delete the next line to include them again.
files <- files[!grepl("_expansion_", basename(files))]
if (length(files) == 0) stop("No *_route_gi.csv found in ", in_dir, " -- run 5b first.")

summary_rows <- list()
persist_rows <- list()

for (f in files) {
  m <- str_match(basename(f), "^(.*?)_(base|anthro)_(rcp45|rcp85)_(.*)_route_gi\\.csv$")
  combo <- tibble(run_label = m[2], model_tag = m[3], scenario = m[4], section = m[5])
  gi <- read.csv(f)
  cat("\n", basename(f), "\n", sep = "")

  for (attr in c("rate", "magnitude")) {
    if (attr == "rate") {
      d <- gi %>% filter(n_category >= min_n_category); vcol <- "pct_mismatch"
    } else {
      d <- gi %>% filter(n_mismatch >= 1); vcol <- "mean_trend_mismatch"
    }
    if (nrow(d) <= max(k_values) + 1) next

    zmat <- matrix(NA_real_, nrow(d), length(k_values), dimnames = list(NULL, paste0("k", k_values)))
    qmat <- zmat
    for (j in seq_along(k_values)) {
      res <- gi_star_at_k(d, vcol, k_values[j])
      zmat[, j] <- res$z
      qmat[, j] <- res$q
      east <- d$longitude > -100
      summary_rows[[length(summary_rows) + 1]] <- bind_cols(
        combo, tibble(attribute = attr, k = k_values[j], n_routes = nrow(d),
          n_hot_95 = sum(res$z >= z_hot), n_hot_99 = sum(res$z >= z_hot99),
          n_hot_fdr = sum(res$q < fdr_q & res$z > 0),
          n_cold_95 = sum(res$z <= -z_hot),
          median_kth_km_east = median(res$kth_km[east]),
          median_kth_km_west = median(res$kth_km[!east])))
    }

    hot  <- zmat >=  z_hot
    cold <- zmat <= -z_hot
    hot_fdr <- qmat < fdr_q & zmat > 0
    ref  <- paste0("k", k_reference)
    persist_rows[[length(persist_rows) + 1]] <- bind_cols(
      combo, tibble(attribute = attr, route = d$route, latitude = d$latitude,
                    longitude = d$longitude,
                    n_category = d$n_category, n_mismatch = d$n_mismatch,
                    pct_mismatch = d$pct_mismatch,
                    mean_trend_mismatch = d$mean_trend_mismatch,
                    n_groups_mismatch = d$n_groups_mismatch,
                    n_k_hot = rowSums(hot), n_k_cold = rowSums(cold),
                    hot_at_ref = hot[, ref],
                    n_k_fdr_hot = rowSums(hot_fdr),
                    fdr_hot_at_ref = hot_fdr[, ref]))
  }
}

summary_df <- bind_rows(summary_rows)
persist_df <- bind_rows(persist_rows)

# Retention: of the reference-k hot spots, what share stay hot at each k, and
# how many routes are hot at EVERY k (the scale-robust core).
ref_stats <- persist_df %>%
  group_by(run_label, model_tag, scenario, section, attribute) %>%
  summarise(n_hot_ref = sum(hot_at_ref),
            n_fdr_hot_ref = sum(fdr_hot_at_ref),
            n_fdr_hot_all_k = sum(n_k_fdr_hot == length(k_values)),
            n_hot_all_k = sum(n_k_hot == length(k_values)),
            n_hot_any_k = sum(n_k_hot >= 1),
            n_cold_all_k = sum(n_k_cold == length(k_values)),
            .groups = "drop")
summary_df <- summary_df %>%
  left_join(ref_stats, by = c("run_label", "model_tag", "scenario", "section", "attribute"))

write.csv(summary_df, file.path(out_dir, "gi_k_sensitivity_summary.csv"), row.names = FALSE)
write.csv(persist_df, file.path(out_dir, "gi_k_sensitivity_route_persistence.csv"), row.names = FALSE)
cat("\nWrote summary + persistence CSVs to", out_dir, "\n")

# Final candidate list: routes that are FDR-significant hot spots (q < fdr_q)
# at EVERY tested k -- robust to both the multiple-testing problem and the
# neighborhood-size choice.
robust_df <- persist_df %>%
  filter(n_k_fdr_hot == length(k_values)) %>%
  arrange(run_label, model_tag, scenario, section, attribute, desc(n_mismatch))
write.csv(robust_df, file.path(out_dir, "robust_refugia_candidates.csv"), row.names = FALSE)
cat("Wrote", nrow(robust_df), "robust (FDR-significant at all k) route rows to robust_refugia_candidates.csv\n")

# Persistence maps: color = number of k values (out of length(k_values)) at
# which the route is an FDR-significant hot spot, rate and magnitude side by side.
pal <- c("grey85", colorRampPalette(c("#FEE0D2", "#67000D"))(length(k_values)))
names(pal) <- 0:length(k_values)

combos <- persist_df %>% distinct(run_label, model_tag, scenario, section)
for (i in seq_len(nrow(combos))) {
  cb <- combos[i, ]
  d <- persist_df %>%
    semi_join(cb, by = c("run_label", "model_tag", "scenario", "section")) %>%
    mutate(attribute = factor(attribute, c("rate", "magnitude"),
                              c("Mismatch rate", "Refugia strength (mean trend)")),
           n_k_fdr_hot = factor(n_k_fdr_hot, 0:length(k_values))) %>%
    arrange(n_k_fdr_hot)
  p <- ggplot() +
    geom_sf(data = us_states, fill = "grey97", color = "grey70", linewidth = 0.25) +
    geom_point(data = d, aes(longitude, latitude, color = n_k_fdr_hot), size = 1.3, alpha = 0.9) +
    scale_color_manual(values = pal, drop = FALSE,
                       name = paste0("# of k values\nFDR-hot (of ", length(k_values), ")")) +
    facet_wrap(~attribute) +
    coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
    labs(title = paste0("Gi* hot-spot persistence across k = ", paste(k_values, collapse = ", "),
                        " (", cb$model_tag, ", ", toupper(sub("rcp", "RCP ", cb$scenario)), ")"),
         subtitle = cb$section, x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(size = 13, face = "bold"))
  fn <- file.path(plot_dir, paste0(cb$run_label, "_", cb$model_tag, "_", cb$scenario, "_",
                                   cb$section, "_k_persistence_map.png"))
  ggsave(fn, p, width = 14, height = 6, dpi = 150, bg = "white")
}
cat("Saved", nrow(combos), "persistence maps to", plot_dir, "\n")
