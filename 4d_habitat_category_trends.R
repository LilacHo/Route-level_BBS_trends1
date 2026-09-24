# 4d_habitat_category_trends.R
#
# Tests the core hypothesis this project is built around: route-level BBS
# population trend should be INCREASING in areas the SDM classifies as
# habitat EXPANSION and DECREASING in areas classified as CONTRACTION
# (Stable falling in between). Answers two questions:
#
#   1. Does the trend follow that pattern -- overall, and separately within
#      each of the 12 bird groups in data/spp_names_codes_group_aou.csv?
#   2. Does the answer change between the "base" and "anthro" SDM/trend
#      models?
#
# Two ways of averaging are computed side by side, because they can and do
# disagree (see the project record for the full worked example):
#   - unweighted_mean  -- simple mean of each SPECIES' own mean trend, every
#                         species counted once regardless of sample size.
#                         Sensitive to species estimated from only 1-2
#                         routes (a single route's idiosyncratic trend can
#                         swing that species' "mean" by +/-15%/yr), which
#                         then counts exactly as much as a 300-route species.
#   - weighted_mean    -- each species' mean trend weighted by how many
#                         routes it was estimated from (Contraction/Stable/
#                         Expansion n). Approximates the average trend per
#                         ROUTE rather than per species, and is far more
#                         stable to the inclusion/exclusion of low-n species
#                         (checked directly: filtering to n > 0/3/5/10 barely
#                         moves the weighted mean but moves the unweighted
#                         mean by several tenths of a %/yr).
# Neither is "more correct" in the abstract -- unweighted treats every
# species' fate as equally important regardless of range size; weighted
# answers "what does a typical monitored ROUTE look like." Given how much
# the unweighted mean moves with an arbitrary n-cutoff, weighted is the more
# robust number to lead with; unweighted is kept alongside as the
# species-equal-weight view.
#
# Reads:  output/files/all_species_habitat_contraction_stable_expansion_2010_2025_base.xlsx
#         output/files/all_species_habitat_contraction_stable_expansion_2010_2025_anthro.xlsx
#         Each workbook has one sheet per scenario (RCP4.5, RCP8.5) and one
#         row per species: Contraction/Stable/Expansion n (route count),
#         mean, sd, and sig (a compact-letter-display code from a per-species
#         pairwise test across the three categories -- shared letters =
#         not significantly different). These are NOT produced by any
#         script in this project (same situation as the hand-built
#         gamma1_ci_comparison_all_species_*.csv files documented in
#         helper/group_stats.R's header) -- they summarize the per-route,
#         per-species trend x SDM-category data that
#         3c_add_SDM_covariates.R writes to
#         output/species_routes_covariates/per_species_sdm/*_route_trends_sdm.csv
#         and 4c_statistical_analysis_and_visualization_covariates.R's PART 4
#         already tests species-by-species. If these workbooks are ever
#         regenerated directly by a script, point in_files below at that
#         script's output instead.
#         data/spp_names_codes_group_aou.csv -- species -> bird Group.
#
# A fourth pseudo-category, "Overall", is added alongside Contraction/
# Stable/Expansion: each species' own overall mean trend across ALL its
# categorized routes, recovered by taking a route-count-weighted average of
# that species' (up to three) category means -- since the categories are a
# non-overlapping partition of a species' routes, this exactly reconstructs
# the grand per-species mean (no need to go back to the route-level data).
# Averaging that Overall value across species the same way as the other
# three categories (unweighted / route-count-weighted) gives each group a
# baseline trend to compare its Contraction/Stable/Expansion values against.
#
# Writes: output/files/habitat_category_group_trend_summary_2010_2025.csv
#         Long format: one row per model x scenario x Group (plus a pooled
#         "All species" row) x category (Contraction/Stable/Expansion/
#         Overall): n_species, total_routes, unweighted_mean, weighted_mean.
#         output/files/habitat_category_group_trend_comparison_2010_2025.csv
#         Wide format (one row per model x scenario x Group): the four
#         category means side by side under both averaging schemes, plus
#         whether Expansion > Contraction under each -- this is the table
#         printed to console below, saved so it doesn't need retyping.

library(here)
library(tidyverse)
library(readxl)

here::i_am("4d_habitat_category_trends.R")

# Settings -------------------------------------------------------------------
in_files <- c(
  base   = here::here("output", "files", "all_species_habitat_contraction_stable_expansion_2010_2025_base.xlsx"),
  anthro = here::here("output", "files", "all_species_habitat_contraction_stable_expansion_2010_2025_anthro.xlsx")
)
group_lookup_csv <- here::here("data", "spp_names_codes_group_aou.csv")
out_csv <- here::here("output", "files", "habitat_category_group_trend_summary_2010_2025.csv")
out_csv_wide <- here::here("output", "files", "habitat_category_group_trend_comparison_2010_2025.csv")

missing_files <- in_files[!file.exists(in_files)]
if (length(missing_files) > 0) {
  stop("Missing input workbook(s): ", paste(missing_files, collapse = ", "),
       " -- see this script's header comment for where these come from.")
}

# Read + combine both workbooks (both RCP sheets each) -----------------------
read_one_workbook <- function(path, model_tag) {
  sheets <- excel_sheets(path)
  map_dfr(sheets, function(sh) {
    d <- read_excel(path, sheet = sh)
    d$scenario <- sh
    d$model <- model_tag
    d
  })
}

all_raw <- map2_dfr(in_files, names(in_files), read_one_workbook)
names(all_raw) <- make.names(names(all_raw))  # "Contraction n" -> "Contraction.n" etc.

# Bird group lookup, de-duplicated the same way every other script in this
# project resolves species that are dual-listed under two groups (e.g.
# Willet under both "coastal" and "waterbirds") -- keep the first-listed
# group only, same convention as 1c/2c/helper/gamma_lookup.R etc.
spp_group <- read.csv(group_lookup_csv, stringsAsFactors = FALSE) %>%
  distinct(Common.Name, Code, .keep_all = TRUE) %>%
  select(Species = Common.Name, Group)

long <- all_raw %>%
  select(Species, model, scenario,
         Contraction_n = Contraction.n, Contraction_mean = Contraction.mean,
         Stable_n = Stable.n, Stable_mean = Stable.mean,
         Expansion_n = Expansion.n, Expansion_mean = Expansion.mean) %>%
  left_join(spp_group, by = "Species")

n_unmatched <- long %>% filter(is.na(Group)) %>% distinct(Species) %>% nrow()
if (n_unmatched > 0) {
  message(n_unmatched, " species had no Group match in ", group_lookup_csv,
          " -- excluded from the by-group summary (still included in All species).")
}

# Each species' own overall mean trend across every categorized route, and
# the route count backing it -- a route-count-weighted average of that
# species' (up to 3) category means exactly recovers the grand per-species
# mean, since Contraction/Stable/Expansion routes are non-overlapping.
long <- long %>%
  rowwise() %>%
  mutate(
    Overall_n = sum(c(Contraction_n, Stable_n, Expansion_n), na.rm = TRUE),
    Overall_mean = {
      n_vec <- c(Contraction_n, Stable_n, Expansion_n)
      m_vec <- c(Contraction_mean, Stable_mean, Expansion_mean)
      keep  <- !is.na(n_vec) & !is.na(m_vec) & n_vec > 0
      if (any(keep)) sum(m_vec[keep] * n_vec[keep]) / sum(n_vec[keep]) else NA_real_
    }
  ) %>%
  ungroup()

# Reshape to one row per model x scenario x species x category ---------------
long_cat <- long %>%
  pivot_longer(cols = c(Contraction_mean, Stable_mean, Expansion_mean, Overall_mean),
               names_to = "category", values_to = "trend_mean") %>%
  mutate(category = sub("_mean$", "", category),
         n_col = case_when(category == "Contraction" ~ Contraction_n,
                           category == "Stable"      ~ Stable_n,
                           category == "Expansion"   ~ Expansion_n,
                           category == "Overall"     ~ Overall_n)) %>%
  filter(!is.na(trend_mean)) %>%
  mutate(category = factor(category, levels = c("Contraction", "Stable", "Expansion", "Overall")))

summarise_trends <- function(df) {
  df %>%
    summarise(n_species = n(), total_routes = sum(n_col, na.rm = TRUE),
              unweighted_mean = mean(trend_mean, na.rm = TRUE),
              weighted_mean   = weighted.mean(trend_mean, n_col, na.rm = TRUE),
              .groups = "drop")
}

by_group <- long_cat %>%
  filter(!is.na(Group)) %>%
  group_by(model, scenario, Group, category) %>%
  summarise_trends()

all_species <- long_cat %>%
  group_by(model, scenario, category) %>%
  summarise_trends() %>%
  mutate(Group = "All species", .after = scenario)

category_summary <- bind_rows(all_species, by_group) %>%
  mutate(Group = factor(Group, levels = c("All species", sort(unique(by_group$Group))))) %>%
  arrange(model, scenario, Group, category)

if (!dir.exists(dirname(out_csv))) dir.create(dirname(out_csv), recursive = TRUE)
write.csv(category_summary, out_csv, row.names = FALSE)
cat("Wrote:", out_csv, "\n")

# Wide comparison table: does Expansion trend beat Contraction trend, per
# model x scenario x Group ("All species" included as its own row), with
# each group's Overall (unw./wtd.) trend alongside as a baseline? -----------
wide_comparison <- category_summary %>%
  pivot_wider(names_from = category,
              values_from = c(n_species, total_routes, unweighted_mean, weighted_mean),
              names_glue = "{category}_{.value}") %>%
  mutate(
    n_species_total = Overall_n_species,
    unweighted_follows_hypothesis = Expansion_unweighted_mean > Contraction_unweighted_mean,
    weighted_follows_hypothesis   = Expansion_weighted_mean   > Contraction_weighted_mean
  ) %>%
  arrange(model, scenario, desc(Group == "All species"), desc(n_species_total))

# Same wide layout, written to its own CSV (one sheet's worth of rows per
# model x scenario, stacked with model/scenario as leading columns) so the
# console table above doesn't have to be re-typed by hand to share/reuse it.
wide_comparison_out <- wide_comparison %>%
  select(model, scenario, Group, n_species = n_species_total,
         `Contraction (unw.)` = Contraction_unweighted_mean,
         `Stable (unw.)`      = Stable_unweighted_mean,
         `Expansion (unw.)`   = Expansion_unweighted_mean,
         `Overall (unw.)`     = Overall_unweighted_mean,
         `Hyp. holds? (unw.)` = unweighted_follows_hypothesis,
         `Contraction (wtd.)` = Contraction_weighted_mean,
         `Stable (wtd.)`      = Stable_weighted_mean,
         `Expansion (wtd.)`   = Expansion_weighted_mean,
         `Overall (wtd.)`     = Overall_weighted_mean,
         `Hyp. holds? (wtd.)` = weighted_follows_hypothesis)

write.csv(wide_comparison_out, out_csv_wide, row.names = FALSE)
cat("Wrote:", out_csv_wide, "\n")

cat("\n=== Per model x scenario x group: category trends vs. that group's own Overall trend ===\n")
for (m in unique(wide_comparison$model)) {
  for (sc in unique(wide_comparison$scenario)) {
    cat("\n---", m, "|", sc, "---\n")
    print(as.data.frame(
      wide_comparison %>%
        filter(model == m, scenario == sc) %>%
        select(Group, n_species_total,
               Contraction_unweighted_mean, Stable_unweighted_mean, Expansion_unweighted_mean,
               Overall_unweighted_mean, unweighted_follows_hypothesis,
               Contraction_weighted_mean, Stable_weighted_mean, Expansion_weighted_mean,
               Overall_weighted_mean, weighted_follows_hypothesis)
    ), digits = 3, row.names = FALSE)
  }
}

cat("\nDone.\n")
