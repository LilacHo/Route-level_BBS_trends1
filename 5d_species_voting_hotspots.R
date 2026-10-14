# 5d_species_voting_hotspots.R
#
# A third way of finding refugia hotspots, on top of 5_hotspots.R (pooled
# mismatch rate/magnitude) and 5b/5c (Getis-Ord Gi* on those pooled values).
# Both of those pool every species together into one number per route, which
# 5b's follow-up check showed has a real failure mode: a single widespread
# species (Eurasian Collared-Dove) could single-handedly drive most of a
# "hot spot" cluster, because pooling lets one species contribute to the
# rate/magnitude at dozens of routes at once.
#
# This script instead runs Gi* SEPARATELY within each eligible species' own
# set of contraction-predicted routes, and lets each qualifying species cast
# AT MOST ONE VOTE, for its single strongest ("least decline") route. A
# route's final score is how many independent species voted for it. Because
# no species can cast more than one vote no matter how many routes it has,
# a widespread species can no longer dominate the result the way it did in
# the pooled analysis -- a high score requires several DIFFERENT species to
# each independently pick the same route as their own best refugium
# candidate, which is a much harder coincidence to produce with one
# ubiquitous species or a run of sparse detections.
#
# Species eligibility (both bounds informed by this project's own data --
# see the console output of this script's first run for the actual
# distribution):
#   - Too FEW contraction-predicted routes and there's not enough spatial
#     information for a species' own Gi* to mean anything (k nearest
#     neighbors needs a real neighborhood to draw from). min_routes_per_species
#     enforces this floor.
#   - Too MANY and the species is (a) common/widespread enough that refugia
#     identification is less conservation-relevant for it specifically, and
#     (b) structurally more likely to rack up a vote just from having more
#     routes to test (Gi* flags roughly a fixed PROPORTION of routes, so
#     more routes means more raw chances to clear the bar) -- neither reason
#     is about data quality, both are about not letting a handful of very
#     common generalists dominate the vote tally. max_routes_per_species
#     enforces this ceiling. In this project's base/RCP4.5 data, a floor of
#     20 and ceiling of 150 keeps 163 of 424 contraction-bearing species
#     eligible, excluding the 60 most ubiquitous (Mallard, Song Sparrow,
#     Canada Goose, ...) and the ~200 with too few routes to test.
#
# Ranks on trend_lci (the credible lower bound), not raw trend -- see the
# project record/methods note: a route's own Gi* result should reflect
# whether the model is actually confident the species isn't declining
# there, not just a point estimate that a single lucky detection can move a
# long way (5b's raw-count check found several "hot" routes driven by a
# single nonzero year in over a decade of surveys).
#
# A voting species' one vote requires that its own within-species Gi* find
# an FDR-significant (q < fdr_q) hot route -- a species with no significant
# local structure in its own contraction range casts no vote at all, rather
# than being forced to vote for whichever route happens to rank highest.
#
# Two scopes, both run here so they can be compared rather than choosing one
# up front (see the project record for why pooling every group together can
# miss the strongest evidence -- cross-guild agreement -- while a single
# group can't detect it at all):
#   - "all_groups" -- every eligible species pooled, so a route can only
#     score high if DIFFERENT species (possibly from different habitat
#     groups) independently agree on it -- the strongest, rarer signal.
#   - one run per Group (aridlands, boreal_forests, ... from
#     data/spp_names_codes_group_aou.csv) -- agreement within one guild,
#     which is more attributable to a shared, nameable mechanism (the same
#     habitat feature) but can't show cross-guild convergence.
#
# A SEPARATE question from "does any one route beat chance": does a
# NEIGHBORHOOD of routes collectively attract more votes than expected, even
# if no single route within it individually clears the permutation line?
# The per-route permutation test above only asks about one route's own
# count in isolation. To ask about neighborhoods, every candidate route
# (every route ANY eligible species could have voted for, not just the ones
# that received a vote) is assigned its score (0 if it got no votes), and
# Getis-Ord Gi* -- the same statistic 5b/5c use, applied here to the vote
# score instead of a pooled mismatch rate -- looks for routes whose k
# nearest neighbors collectively carry more votes than chance. A cluster of
# several adjacent routes with 1 vote each can show up here even though none
# of them individually stood out above.
#
# Score significance (per-route): a permutation null. For each of n_permutations
# iterations, every ACTUAL voting species is reassigned a uniformly random
# route from its own eligible contraction-route set (not just its hot
# routes) instead of its real pick, and the resulting scores are tallied.
# This asks: if these same species' "best" routes had been randomly located
# within their own ranges instead of the model's actual estimates, how often
# would that many species land on the same route by chance? The 95th
# percentile of the permuted max score is reported as score_threshold_p05 --
# an observed score above that line is a meaningfully unlikely coincidence,
# not just "the most votes in this particular run."
#
# Reads:  output/species_routes_covariates/per_species_sdm/*_route_trends_sdm.csv
#         (written by 3c_add_SDM_covariates.R) -- same input as 5_hotspots.R/5b
# Writes: output/species_routes_covariates/hotspots/species_voting/<model_tag>_<scenario>_<route_scope>_<rank_col>_<voting_rule>_<scope>_mindetect<N>_votes.csv
#         output/species_routes_covariates/hotspots/species_voting/<model_tag>_<scenario>_<route_scope>_<rank_col>_<voting_rule>_<scope>_mindetect<N>_route_scores.csv
#         output/species_routes_covariates/hotspots/species_voting/<model_tag>_<scenario>_<route_scope>_<rank_col>_<voting_rule>_<scope>_mindetect<N>_neighborhood_gi.csv
#         output/species_routes_covariates/hotspots/species_voting/<model_tag>_<scenario>_summary_across_rules_and_scopes_mindetect<N>.csv
#         output/species_routes_covariates/hotspots/species_voting/plots/<model_tag>_<scenario>_<route_scope>_<rank_col>_<voting_rule>_<scope>_mindetect<N>_vote_map.png
#         output/species_routes_covariates/hotspots/species_voting/plots/<model_tag>_<scenario>_<route_scope>_<rank_col>_<voting_rule>_<scope>_mindetect<N>_neighborhood_map.png

library(here)
library(tidyverse)
library(sf)
library(spdep)

here::i_am("5d_species_voting_hotspots.R")

# Settings ---------------------------------------------------------------
# Kept to one combination by default so a first run stays reviewable; extend
# to c("base","anthro") / c("rcp45","rcp85") once the approach is settled.
model_tags <- c("anthro")
scenarios  <- c("rcp45")

# Compare all vote-weighting rules in one run -- see the header comment
# above species_votes() for what each does. Narrowed to "top1" for the
# route_scopes comparison below (each rule x route_scope combination is a
# full separate run); re-add "top3"/"fractional"/"cluster" once a route_scope
# is settled on.
voting_rules <- c("top1", "top3", "fractional", "cluster")

# Which per-route value each species' own Gi* is run on: "trend_lci" (the
# credible lower bound -- the project's primary choice, penalizes uncertain
# estimates) or "trend" (the posterior mean point estimate, ignores
# uncertainty). Both are run so the choice's effect can be compared.
rank_cols <- c("trend_lci")

# Minimum number of years (of first_year:last_year, on standard run_type == 1
# route-years) in which the species was actually DETECTED at a route for that
# route to be eligible to receive its vote. A "credible" trend can be built
# from almost no local data -- the iCAR model borrows strength from
# neighboring routes -- and a raw-count check of this project's hot spots
# found many votes resting on 1-2 detection years (ranking on trend_lci
# alone still left ~41% of cluster votes on <=2 detection years). Applied at
# the VOTE level: a species' hot routes are first dropped if they fall below
# this, and the top1/top3/cluster rules then choose among what remains; the
# permutation null's placement pools are restricted the same way so the test
# stays like-for-like. 0 disables the filter (and skips loading the raw data).
min_detect_years <- 0   # disabled (0) -- did not prove useful; see report
first_year <- 2010
last_year  <- 2025

# Only let species vote if their own SDM has demonstrably worked: in
# output/.../per_species_sdm_stats/species_effect_sizes_<model>_<RCP>.csv
# (the contraction-vs-expansion trend comparison), the species must be
# FDR-significant (sig_global_fdr) AND in the SDM-consistent direction
# (r_rb > 0, i.e. expansion-predicted routes trend higher than contraction-
# predicted ones -- the sign convention that matches the pooled result).
# Rationale: for a species whose SDM shows no detectable link to observed
# trends (about 3/4 of species), a route "beating the prediction" is just
# noise around no relationship, so it says nothing about refugia.
require_sdm_skill <- FALSE   # disabled -- did not prove useful; see report

# "contraction_only" (the original design): each species' own Gi* runs only
# on routes the SDM predicts CONTRACTION for -- a hotspot there specifically
# means "defying a decline prediction." "all_routes" runs the identical
# procedure on EVERY route the species has (any SDM category), regardless of
# what's predicted there -- a validity check: if the same routes come out on
# top either way, they may just be generically good/resilient routes for
# that species, not specifically refugia from a predicted contraction. If
# different routes come out, that's evidence the contraction restriction is
# doing real, meaningful work rather than just relabeling ordinary variation.
route_scopes <- c("contraction_only")

require_route_converged <- TRUE

min_routes_per_species <- 20    # floor -- enough routes for k neighbors to mean something
max_routes_per_species <- 300   # ceiling -- excludes the handful of most-ubiquitous species.
# Raised from 150: at floor=20, the eligible-species count vs. ceiling curve
# flattens hard around 300-400 (163 eligible at 150 -> 202 at 300 -> 223 with
# no ceiling at all), so 150 sat at the 85th percentile of route-count and
# excluded a long thin tail of genuinely-common-but-not-ubiquitous species
# (Bobolink, Wood Thrush, Bald Eagle, ...) without excluding many MORE than
# 300 does. The Mallard/Song Sparrow/Canada Goose tier (400+ routes) this
# ceiling exists to exclude stays excluded either way.

k_neighbors <- 8      # same convention as 5b/5c, applied within each species' OWN route set
fdr_q       <- 0.05   # Benjamini-Hochberg threshold for a species' own Gi* (same as 5b)

n_permutations <- 499
set.seed(20260101)    # reproducible permutation null

in_dir  <- here::here("output", "species_routes_covariates", "per_species_sdm")
out_dir <- here::here("output", "species_routes_covariates", "hotspots", "species_voting")
plot_dir <- file.path(out_dir, "plots")
if (!dir.exists(out_dir))  dir.create(out_dir,  recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Read every SDM CSV (same input as 5_hotspots.R / 5b) ----------------------
sdm_files <- list.files(in_dir, pattern = "_route_trends_sdm\\.csv$", full.names = TRUE)
if (length(sdm_files) == 0) {
  stop("No SDM CSVs found in ", in_dir, " — run 3c_add_SDM_covariates.R first.")
}
cat("\nCombining all species SDM files...\n")
all_sdm_raw <- sdm_files %>% map_dfr(read.csv)

if (require_route_converged) {
  n_before <- nrow(all_sdm_raw)
  all_sdm_raw <- all_sdm_raw %>% filter(route_converged == TRUE)
  cat("require_route_converged = TRUE -> dropped", n_before - nrow(all_sdm_raw),
      "of", n_before, "row(s).\n")
}

group_category <- function(x) {
  dplyr::case_when(
    x %in% c(1, 2, 3) ~ "Contraction",
    x == 4            ~ "Stable",
    x %in% c(5, 6, 7) ~ "Expansion",
    TRUE              ~ NA_character_
  )
}

# Detection-year lookup: for every (state_num, route number, AOU code), in how
# many of first_year:last_year the species was detected (species_total > 0) on
# a standard (run_type == 1) run, from the raw BBS data bbsBayes2 already
# caches locally. species_code -> AOU comes from this project's own
# data/spp_names_codes_group_aou.csv (spp.num), so no name matching is needed.
if (min_detect_years > 0) {
  cat("\nLoading cached raw BBS data to count detection years per species-route...\n")
  raw_bbs <- bbsBayes2::load_bbs_data(level = "state")
  good_routes <- raw_bbs$routes %>%
    filter(run_type == 1, year >= first_year, year <= last_year) %>%
    distinct(route_data_id, state_num, route, year)
  detect_lookup <- raw_bbs$birds %>%
    semi_join(good_routes, by = c("route_data_id", "state_num", "route", "year")) %>%
    filter(species_total > 0) %>%
    group_by(state_num, route, aou) %>%
    summarise(n_detect_years = n_distinct(year), .groups = "drop") %>%
    rename(route_num = route)
  rm(raw_bbs, good_routes)
  species_aou <- read.csv(here::here("data", "spp_names_codes_group_aou.csv")) %>%
    distinct(species_code = Code, aou = spp.num)
}

data("us_states", package = "spData")
us_states <- st_transform(us_states, 4326)
us_bbox <- st_bbox(us_states)
map_xlim <- c(us_bbox["xmin"] - 1, us_bbox["xmax"] + 1)
map_ylim <- c(us_bbox["ymin"] - 1, us_bbox["ymax"] + 1)

# Haversine distance in km -- used only by the "cluster" voting rule to
# decide whether two of a species' own hot routes represent the same
# candidate refugium or two separate ones.
hav_km <- function(lon1, lat1, lon2, lat2) {
  r <- 6371; p <- pi / 180
  a <- sin((lat2 - lat1) * p / 2)^2 +
    cos(lat1 * p) * cos(lat2 * p) * sin((lon2 - lon1) * p / 2)^2
  2 * r * asin(pmin(1, sqrt(a)))
}

# Single-linkage clustering of a small set of points by a distance cutoff
# (connected components of the "within cutoff_km" graph) -- no extra package
# needed for the handful of points any one species' hot routes involve.
cluster_by_distance <- function(lon, lat, cutoff_km) {
  n <- length(lon)
  cluster_id <- seq_len(n)
  for (i in seq_len(n - 1)) for (j in (i + 1):n) {
    if (hav_km(lon[i], lat[i], lon[j], lat[j]) <= cutoff_km) {
      old <- cluster_id[j]; new <- cluster_id[i]
      cluster_id[cluster_id == old] <- new
    }
  }
  match(cluster_id, unique(cluster_id))
}

# Distance within which two of a species' own hot routes count as the same
# candidate refugium under voting_rule == "cluster" -- roughly the k=8
# neighbor scale found elsewhere in this project (40-100 km), with margin.
cluster_distance_km <- 150

# One species' own Gi* on trend_lci within its own contraction routes, then
# turns its FDR-significant hot routes into 0, 1, or several votes depending
# on voting_rule:
#   "top1"       -- one vote (weight 1) for its single strongest hot route.
#   "top3"       -- up to 3 full-weight (1 each) votes, for its 3 strongest
#                   hot routes by z (fewer if it has fewer than 3). Unlike
#                   "cluster" below, this caps the COUNT of votes but not
#                   their spatial spread -- 3 routes that are actually part
#                   of the same local patch still count as 3 separate votes.
#   "fractional" -- a vote for EVERY hot route, each weighted 1/(# hot routes)
#                   so the species' total contribution is still capped at 1.
#   "cluster"    -- hot routes within cluster_distance_km of each other are
#                   treated as one candidate refugium; one full-weight (1)
#                   vote per resulting cluster, at that cluster's strongest
#                   route. A species can cast >1 vote here, but only for
#                   locations far enough apart to be independent evidence.
# Returns NULL if the species doesn't clear min_routes_per_species, or casts
# no vote (no FDR-significant hot route at all).
species_votes <- function(sp_df, voting_rule, rank_col) {
  if (nrow(sp_df) < max(k_neighbors + 1, min_routes_per_species)) return(NULL)
  coords <- as.matrix(sp_df[, c("longitude", "latitude")])
  nb <- suppressWarnings(knn2nb(knearneigh(coords, k = k_neighbors)))
  listw <- nb2listw(include.self(nb), style = "B")
  z <- as.numeric(localG(sp_df[[rank_col]], listw))
  q <- p.adjust(2 * pnorm(-abs(z)), method = "BH")
  hot <- which(q < fdr_q & z > 0)
  # Drop hot routes with too few detection years BEFORE any rule picks among
  # them (see min_detect_years above).
  if (min_detect_years > 0) hot <- hot[sp_df$n_detect_years[hot] >= min_detect_years]
  if (length(hot) == 0) return(NULL)

  pick <- function(idx, w) tibble(route = sp_df$route[idx], latitude = sp_df$latitude[idx],
                                  longitude = sp_df$longitude[idx], z = z[idx], q = q[idx],
                                  rank_value_at_route = sp_df[[rank_col]][idx],
                                  n_detect_years = sp_df$n_detect_years[idx], weight = w)

  if (voting_rule == "top1") {
    best <- hot[which.max(z[hot])]
    pick(best, 1)
  } else if (voting_rule == "top3") {
    top_n <- hot[order(z[hot], decreasing = TRUE)][seq_len(min(3, length(hot)))]
    pick(top_n, 1)
  } else if (voting_rule == "fractional") {
    pick(hot, 1 / length(hot))
  } else if (voting_rule == "cluster") {
    cl <- cluster_by_distance(sp_df$longitude[hot], sp_df$latitude[hot], cluster_distance_km)
    reps <- vapply(unique(cl), function(g) hot[cl == g][which.max(z[hot[cl == g]])], integer(1))
    pick(reps, 1)
  } else {
    stop("Unknown voting_rule: ", voting_rule)
  }
}

# Gi* on the vote score across every CANDIDATE route (not just the ones that
# got a vote) -- zero-filled, so a neighborhood of routes with modest scores
# can stand out even when no single route in it beats the per-route
# permutation line above.
neighborhood_score_gi <- function(candidate_routes, route_scores) {
  full <- candidate_routes %>%
    left_join(route_scores %>% select(route, score), by = "route") %>%
    mutate(score = replace_na(score, 0))
  if (nrow(full) < k_neighbors + 1 || sum(full$score > 0) < 2) return(NULL)
  coords <- as.matrix(full[, c("longitude", "latitude")])
  nb <- suppressWarnings(knn2nb(knearneigh(coords, k = k_neighbors)))
  listw <- nb2listw(include.self(nb), style = "B")
  z <- as.numeric(localG(full$score, listw))
  q <- p.adjust(2 * pnorm(-abs(z)), method = "BH")
  full %>% mutate(gi_z = z, gi_q = q, gi_hot = gi_q < fdr_q & gi_z > 0)
}

run_voting <- function(input_df, scope_label, voting_rule, rank_col) {
  route_counts <- input_df %>% count(species_code, name = "n_routes")
  eligible_codes <- route_counts %>%
    filter(n_routes >= min_routes_per_species, n_routes <= max_routes_per_species) %>%
    pull(species_code)

  cat("\n--- Scope:", scope_label, "| voting_rule:", voting_rule, "| ranked on:", rank_col, "---\n")
  cat("  Species with >=1 route in this route_scope:", n_distinct(input_df$species_code), "\n")
  cat("  Eligible (", min_routes_per_species, "-", max_routes_per_species, " routes):",
      length(eligible_codes), "\n", sep = "")

  if (length(eligible_codes) < 2) {
    cat("  Fewer than 2 eligible species -- skipping this scope.\n")
    return(NULL)
  }

  eligible_df <- input_df %>% filter(species_code %in% eligible_codes)

  votes <- list()
  for (code in eligible_codes) {
    sp_df <- eligible_df %>% filter(species_code == code)
    v <- species_votes(sp_df, voting_rule, rank_col)
    if (!is.null(v)) {
      votes[[code]] <- v %>% mutate(species_code = code,
                                    species = sp_df$species[1], group = sp_df$group[1])
    }
  }
  votes_df <- bind_rows(votes)
  cat("  Species that cast >=1 vote (found an FDR-significant hot route):", n_distinct(votes_df$species_code), "\n")
  cat("  Total (species, route) votes cast:", nrow(votes_df), "\n")

  if (nrow(votes_df) == 0) {
    cat("  No species voted -- skipping this scope.\n")
    return(NULL)
  }

  route_scores <- votes_df %>%
    group_by(route) %>%
    summarise(latitude = first(latitude), longitude = first(longitude),
              score = sum(weight), n_species = n_distinct(species_code),
              species_list = paste(species, collapse = "; "),
              groups_list = paste(unique(group), collapse = "; "),
              n_groups = n_distinct(group), .groups = "drop") %>%
    arrange(desc(score))

  # Permutation null: reassign each ACTUAL voting species' vote(s) to
  # uniformly random, DISTINCT route(s) from its OWN eligible contraction-
  # route set (not just its hot routes) -- same number of votes and same
  # weights as its real picks, just relocated -- retally, repeat. Tests
  # whether the observed agreement between species is more than random
  # placement within each species' own range would produce. Works the same
  # way regardless of voting_rule since it just replays each species' own
  # (possibly multiple, possibly fractional) weight structure at random
  # locations.
  voting_codes <- unique(votes_df$species_code)
  species_weights <- setNames(lapply(voting_codes, function(code) votes_df$weight[votes_df$species_code == code]),
                              voting_codes)
  # Placement pool = the species' routes that COULD have received its vote
  # (>= min_detect_years detections), so random placement is held to the same
  # standard as the real picks.
  species_pools <- setNames(lapply(voting_codes, function(code) {
    sp_rows <- eligible_df[eligible_df$species_code == code, ]
    ok <- if (min_detect_years > 0) sp_rows$n_detect_years >= min_detect_years else rep(TRUE, nrow(sp_rows))
    pool <- unique(sp_rows$route[ok])
    if (length(pool) == 0) unique(sp_rows$route) else pool
  }), voting_codes)

  perm_max <- numeric(n_permutations)
  for (p in seq_len(n_permutations)) {
    perm_routes <- character(0); perm_weights <- numeric(0)
    for (code in voting_codes) {
      w <- species_weights[[code]]; pool <- species_pools[[code]]
      chosen <- sample(pool, min(length(w), length(pool)))
      perm_routes <- c(perm_routes, chosen)
      perm_weights <- c(perm_weights, w[seq_along(chosen)])
    }
    perm_max[p] <- max(tapply(perm_weights, perm_routes, sum))
  }
  score_threshold_p05 <- as.numeric(quantile(perm_max, 0.95))
  cat("  Permutation null (n=", n_permutations, "): max score at random placement, ",
      "median = ", round(median(perm_max), 2), ", 95th pct = ", round(score_threshold_p05, 2), "\n", sep = "")
  cat("  Observed top score:", round(max(route_scores$score), 2),
      "| routes at/above the p<0.05 threshold:",
      sum(route_scores$score > score_threshold_p05), "\n")

  # Neighborhood-level check: does a group of adjacent routes collectively
  # carry more votes than expected, even if no single route qualifies above?
  candidate_routes <- eligible_df %>% distinct(route, latitude, longitude)
  nbhd <- neighborhood_score_gi(candidate_routes, route_scores)
  if (!is.null(nbhd)) {
    n_nbhd_hot <- sum(nbhd$gi_hot)
    cat("  Neighborhood Gi* on vote score (", nrow(candidate_routes), " candidate routes): ",
        n_nbhd_hot, " routes in an FDR-significant high-vote neighborhood\n", sep = "")
  } else {
    cat("  Neighborhood Gi* skipped (too few candidate routes or too few votes).\n")
  }

  list(votes = votes_df %>% mutate(scope = scope_label),
       route_scores = route_scores %>% mutate(scope = scope_label,
                                              score_threshold_p05 = score_threshold_p05,
                                              significant = score > score_threshold_p05),
       neighborhood = if (!is.null(nbhd)) nbhd %>% mutate(scope = scope_label) else NULL,
       n_eligible = length(eligible_codes), n_voting = n_distinct(votes_df$species_code),
       score_threshold_p05 = score_threshold_p05)
}

vote_rule_caption <- c(
  top1       = "Each eligible species casts at most one vote, for its single strongest (FDR-significant) route.",
  top3       = "Each eligible species casts up to 3 full-weight votes, for its 3 strongest FDR-significant routes.",
  fractional = "Each eligible species casts a vote for EVERY FDR-significant route, each weighted 1/(# such routes).",
  cluster    = paste0("Each eligible species casts one full vote per spatially distinct cluster (>", cluster_distance_km,
                      "km apart) of its FDR-significant routes.")
)

make_vote_map <- function(route_scores, scope_label, model_tag, scenario, voting_rule, rank_col, file_path) {
  thr <- route_scores$score_threshold_p05[1]
  p <- ggplot() +
    geom_sf(data = us_states, fill = "grey96", color = "grey65", linewidth = 0.3) +
    geom_point(data = route_scores %>% filter(!significant),
               aes(x = longitude, y = latitude, size = score), color = "grey50", alpha = 0.5) +
    geom_point(data = route_scores %>% filter(significant),
               aes(x = longitude, y = latitude, size = score), color = "#B2182B", alpha = 0.9) +
    scale_size_continuous(name = "Vote\nscore", range = c(1.5, 7)) +
    coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
    labs(title = paste0("Species-voting hotspots (", voting_rule, ", ", rank_col, "): ", scope_label, " (", model_tag, ", ",
                        toupper(sub("rcp", "RCP ", scenario)), ")"),
         subtitle = paste0(nrow(route_scores), " routes with >=1 vote | red = score > ", round(thr, 2),
                           " (permutation p < 0.05)"),
         caption = vote_rule_caption[[voting_rule]],
         x = NULL, y = NULL) +
    theme_minimal() +
    theme(plot.title    = element_text(size = 15, face = "bold"),
          plot.subtitle = element_text(size = 10),
          plot.caption  = element_text(size = 8, color = "grey40"),
          axis.text     = element_text(size = 8))
  ggsave(file_path, p, width = 9, height = 6.5, dpi = 150, bg = "white")
  cat("  Saved plot:", basename(file_path), "\n")
}

make_neighborhood_map <- function(nbhd, scope_label, model_tag, scenario, voting_rule, rank_col, file_path) {
  p <- ggplot() +
    geom_sf(data = us_states, fill = "grey96", color = "grey65", linewidth = 0.3) +
    geom_point(data = nbhd %>% filter(!gi_hot, score == 0), aes(x = longitude, y = latitude),
               color = "grey85", size = 0.6, alpha = 0.5) +
    geom_point(data = nbhd %>% filter(!gi_hot, score > 0), aes(x = longitude, y = latitude, size = score),
               color = "grey45", alpha = 0.6) +
    geom_point(data = nbhd %>% filter(gi_hot), aes(x = longitude, y = latitude, size = score),
               color = "#B2182B", alpha = 0.9) +
    scale_size_continuous(name = "Vote\nscore", range = c(1, 6)) +
    coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
    labs(title = paste0("Neighborhood vote clustering (", voting_rule, ", ", rank_col, "): ", scope_label, " (", model_tag, ", ",
                        toupper(sub("rcp", "RCP ", scenario)), ")"),
         subtitle = paste0(nrow(nbhd), " candidate routes | red = FDR-significant high-vote neighborhood ",
                           "(Gi* on vote score, k=", k_neighbors, ")"),
         caption = "Tests whether NEARBY routes collectively carry more votes than chance, not just one route alone.",
         x = NULL, y = NULL) +
    theme_minimal() +
    theme(plot.title    = element_text(size = 15, face = "bold"),
          plot.subtitle = element_text(size = 10),
          plot.caption  = element_text(size = 8, color = "grey40"),
          axis.text     = element_text(size = 8))
  ggsave(file_path, p, width = 9, height = 6.5, dpi = 150, bg = "white")
  cat("  Saved plot:", basename(file_path), "\n")
}

# Main loop ----------------------------------------------------------------
for (model_tag in model_tags) {
  target_sdm <- all_sdm_raw %>% filter(model == model_tag)

  for (scenario in scenarios) {
    cat("\n============================================================\n")
    cat("model_tag:", model_tag, "| scenario:", scenario, "\n")
    cat("============================================================\n")

    analysis <- target_sdm %>%
      filter(!is.na(.data[[scenario]]), !is.na(trend_lci), !is.na(trend)) %>%
      transmute(species, species_code, group, route, latitude, longitude, trend_lci, trend,
                category = group_category(.data[[scenario]]))

    if (require_sdm_skill) {
      skill_file <- here::here("output", "species_routes_covariates", "per_species_sdm_stats",
                               paste0("species_effect_sizes_", model_tag, "_", toupper(scenario), ".csv"))
      skill <- read.csv(skill_file)
      skilled_codes <- skill$species_code[as.logical(skill$sig_global_fdr) & skill$r_rb > 0]
      cat("SDM-skill filter: ", length(skilled_codes), " of ", nrow(skill),
          " tested species pass (FDR-significant, SDM-consistent); ",
          n_distinct(analysis$species_code[analysis$species_code %in% skilled_codes]),
          " are present in this analysis.\n", sep = "")
      analysis <- analysis %>% filter(species_code %in% skilled_codes)
    }

    if (min_detect_years > 0) {
      analysis <- analysis %>%
        left_join(species_aou, by = "species_code") %>%
        mutate(state_num = as.integer(sub("-.*", "", route)),
               route_num = as.integer(sub(".*-", "", route))) %>%
        left_join(detect_lookup, by = c("state_num", "route_num", "aou")) %>%
        mutate(n_detect_years = replace_na(n_detect_years, 0L)) %>%
        select(-state_num, -route_num, -aou)
    } else {
      analysis$n_detect_years <- NA_integer_
    }
    contraction_df <- analysis %>% filter(category == "Contraction")

    summary_rows <- list()

    for (route_scope in route_scopes) {
    working_df <- if (route_scope == "contraction_only") contraction_df else analysis
    scopes <- c("all_groups", sort(unique(working_df$group)))
    cat("\n>>> route_scope:", route_scope, "(", nrow(working_df), "species-route rows,",
        n_distinct(working_df$species_code), "species ) <<<\n")

    for (rank_col in rank_cols) {
    for (voting_rule in voting_rules) {
    for (scope in scopes) {
      scope_data <- if (scope == "all_groups") working_df else working_df %>% filter(group == scope)
      result <- run_voting(scope_data, scope, voting_rule, rank_col)
      if (is.null(result)) next

      file_prefix <- paste0(model_tag, "_", scenario, "_", route_scope, "_", rank_col, "_",
                            voting_rule, "_", scope, "_mindetect", min_detect_years,
                            if (require_sdm_skill) "_skill" else "")
      write.csv(result$votes, file.path(out_dir, paste0(file_prefix, "_votes.csv")), row.names = FALSE)
      write.csv(result$route_scores, file.path(out_dir, paste0(file_prefix, "_route_scores.csv")), row.names = FALSE)

      if (any(result$route_scores$score >= 2)) {
        make_vote_map(result$route_scores, scope, model_tag, scenario, voting_rule, rank_col,
                      file.path(plot_dir, paste0(file_prefix, "_vote_map.png")))
      } else {
        cat("  No route reached score >= 2 -- skipping map for this scope.\n")
      }

      n_nbhd_hot <- 0
      if (!is.null(result$neighborhood)) {
        write.csv(result$neighborhood, file.path(out_dir, paste0(file_prefix, "_neighborhood_gi.csv")),
                 row.names = FALSE)
        n_nbhd_hot <- sum(result$neighborhood$gi_hot)
        if (n_nbhd_hot > 0) {
          make_neighborhood_map(result$neighborhood, scope, model_tag, scenario, voting_rule, rank_col,
                                file.path(plot_dir, paste0(file_prefix, "_neighborhood_map.png")))
        }
      }

      summary_rows[[paste(route_scope, rank_col, voting_rule, scope)]] <- tibble(
        model_tag, scenario, route_scope, rank_col, voting_rule, scope, min_detect_years,
        n_eligible_species = result$n_eligible, n_voting_species = result$n_voting,
        n_votes_cast = nrow(result$votes),
        n_routes_with_votes = nrow(result$route_scores),
        max_score = max(result$route_scores$score),
        score_threshold_p05 = result$score_threshold_p05,
        n_routes_significant = sum(result$route_scores$significant),
        n_routes_in_hot_neighborhood = n_nbhd_hot)
    }
    }
    }
    }

    if (length(summary_rows) > 0) {
      summary_df <- bind_rows(summary_rows)
      write.csv(summary_df,
               file.path(out_dir, paste0(model_tag, "_", scenario, "_summary_across_rules_and_scopes_mindetect",
                                         min_detect_years, if (require_sdm_skill) "_skill" else "", ".csv")),
               row.names = FALSE)
      cat("\n--- Summary across route scopes, voting rules, and scopes (", model_tag, "/", scenario, ") ---\n", sep = "")
      print(as.data.frame(summary_df))
    }
  }
}

cat("\n=== Species-voting hotspot analysis complete ===\n")
