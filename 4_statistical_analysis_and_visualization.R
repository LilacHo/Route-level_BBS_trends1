# 4_statistical_analysis_and_visualization.R
#
# Combines the route-level trend statistical analysis with violin plots
# annotated with compact letter displays (CLD).
#
# The CLD letters are derived from the same pairwise Wilcoxon (BH-adjusted)
# tests used in the statistical analysis: groups that share a letter are NOT
# significantly different (alpha = 0.05).
#
# Focus response throughout: trend (annual % change in route-level counts,
# written by 2_generate_route_trend_csvs.R and carried through 3_add_SDM.R).

library(here)
library(tidyverse)
library(multcompView) # multcompView converts pairwise p-values into compact letter displays.

here::i_am("4_statistical_analysis_and_visualization.R")

# Settings -----------------------------------------------------------------
land_cover <- "grasslands"   # target group: species are matched to this group
                             # via data/spp_names_codes_group_aou.csv (the
                             # per-species CSVs from 3_add_SDM.R are named by
                             # species, not by land_cover, so filtering by
                             # group requires this lookup rather than a
                             # filename prefix)

out_dir   <- here::here("output", "species_routes_sdm")
stats_dir <- here::here("output", "species_routes_sdm_stats")
plot_dir  <- here::here("output", "species_routes_sdm_plot")
per_species_dir <- file.path(plot_dir, "per_species")
if (!dir.exists(stats_dir))       dir.create(stats_dir,       recursive = TRUE)
if (!dir.exists(plot_dir))        dir.create(plot_dir,        recursive = TRUE)
if (!dir.exists(per_species_dir)) dir.create(per_species_dir, recursive = TRUE)

# Species lookup table (Code -> Group), used to resolve each species CSV's
# land-cover group since the file name itself carries no land_cover prefix.
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv")) %>%
  distinct(Code, Group)

# Read every SDM CSV written by 3_add_SDM.R, then attach each row's land-cover
# group via species_code so it can be filtered per land_cover below.
sdm_files <- list.files(out_dir, pattern = "_route_trends_sdm\\.csv$",
                        full.names = TRUE)

if (length(sdm_files) == 0) {
  stop("No SDM CSVs found in ", out_dir, " — run 3_add_SDM.R first.")
}

cat("\nCombining all species SDM files...\n")

all_sdm_raw <- sdm_files %>%
  map_dfr(read.csv) %>%
  left_join(spp_df, by = c("species_code" = "Code"))

target_sdm <- all_sdm_raw %>% filter(Group == land_cover)

if (nrow(target_sdm) == 0) {
  stop("No SDM rows matched land_cover = '", land_cover, "' after joining ",
       "species groups from spp_names_codes_group_aou.csv")
}

cat("Total rows:", nrow(target_sdm), "\n")
cat(land_cover, "species:", length(unique(target_sdm$species_code)), "\n")

# Build combined analysis data ---------------------------------------------
analysis_45 <- target_sdm %>%
  filter(!is.na(rcp45)) %>%
  transmute(route, species_code, category = rcp45, trend,
            route_num = as.integer(sub("^\\d+-", "", route)))

analysis_85 <- target_sdm %>%
  filter(!is.na(rcp85)) %>%
  transmute(route, species_code, category = rcp85, trend,
            route_num = as.integer(sub("^\\d+-", "", route)))

# Grouped categories: Contraction = 1,2,3 | Stable = 4 | Expansion = 5,6,7
#   (category 0 = never suitable, excluded)
group_category <- function(x) {
  dplyr::case_when(
    x %in% c(1, 2, 3) ~ "Contraction",
    x == 4            ~ "Stable",
    x %in% c(5, 6, 7) ~ "Expansion",
    TRUE              ~ NA_character_
  )
}

analysis_45_grp <- analysis_45 %>%
  mutate(change_group = factor(group_category(category),
                               levels = c("Contraction", "Stable", "Expansion")))

analysis_85_grp <- analysis_85 %>%
  mutate(change_group = factor(group_category(category),
                               levels = c("Contraction", "Stable", "Expansion")))

target_species_list <- unique(target_sdm$species_code)

cat8_levels <- as.character(0:7)
grp_levels  <- c("Contraction", "Stable", "Expansion")


# ==========================================================================
# Functions ####
# ==========================================================================

# Run Kruskal-Wallis + pairwise Wilcoxon (BH-corrected). Prints a report and
# returns the fitted objects (so the pairwise matrix can feed the CLD).
run_category_tests <- function(df, response, group, scenario_label, response_label) {
  d <- df %>%
    filter(!is.na(.data[[response]]), !is.na(.data[[group]])) %>%
    mutate(grp = factor(.data[[group]]))

  cat("\n--------------------------------------------------\n")
  cat(scenario_label, "|", response_label, "by", group, "\n")
  cat("--------------------------------------------------\n")

  # Need at least 2 groups with data
  if (length(unique(d$grp)) < 2) {
    cat("  Not enough groups for testing.\n")
    return(invisible(NULL))
  }

  # Group sizes and summary
  summ_tab <- d %>%
    group_by(grp) %>%
    summarise(n = n(),
              median = median(.data[[response]]),
              mean = mean(.data[[response]]),
              .groups = "drop")
  cat("\nGroup summary:\n")
  print(as.data.frame(summ_tab), row.names = FALSE)

  # Kruskal-Wallis omnibus test
  kw <- kruskal.test(d[[response]] ~ d$grp)
  cat("\nKruskal-Wallis: chi-sq =", round(kw$statistic, 2),
      ", df =", kw$parameter,
      ", p =", format.pval(kw$p.value, digits = 3), "\n")

  # Pairwise Wilcoxon (BH-corrected)
  pw <- pairwise.wilcox.test(d[[response]], d$grp, p.adjust.method = "BH")
  cat("\nPairwise Wilcoxon (BH-adjusted p-values):\n")
  print(round(pw$p.value, 4))

  invisible(list(kruskal = kw, pairwise = pw, data = d))
}

# Build a full symmetric p-value matrix from a pairwise.wilcox.test result.
# pairwise.wilcox.test returns a lower-triangular matrix; multcompLetters
# needs a named vector of p-values for all pairs.
cld_from_pairwise <- function(pw, alpha = 0.05) {
  if (is.null(pw) || is.null(pw$p.value)) return(NULL)

  pmat <- pw$p.value
  groups <- union(rownames(pmat), colnames(pmat))

  # Assemble named vector "g1-g2" = p for every available pair
  pairs <- character(0)
  pvals <- numeric(0)
  for (r in rownames(pmat)) {
    for (cc in colnames(pmat)) {
      val <- pmat[r, cc]
      if (!is.na(val)) {
        pairs <- c(pairs, paste(r, cc, sep = "-"))
        pvals <- c(pvals, val)
      }
    }
  }
  if (length(pvals) == 0) return(NULL)

  names(pvals) <- pairs
  res <- tryCatch(
    multcompView::multcompLetters(pvals, threshold = alpha),
    error = function(e) NULL
  )
  if (is.null(res)) return(NULL)

  tibble(group = names(res$Letters),
         cld   = unname(res$Letters))
}

# Quietly compute the pairwise Wilcoxon (BH) result for CLD purposes only.
# No printing — used by the plotting layer so it doesn't duplicate the
# statistical report already written to the stats files.
pairwise_quiet <- function(df, response, group) {
  d <- df %>%
    filter(!is.na(.data[[response]]), !is.na(.data[[group]])) %>%
    mutate(grp = factor(.data[[group]]))
  if (length(unique(d$grp)) < 2) return(NULL)
  pairwise.wilcox.test(d[[response]], d$grp, p.adjust.method = "BH")
}

# Compute per-group CLD labels for a single scenario/grouping, returning a
# tibble with group, scenario, the letter, and a y position for the label.
compute_cld_layer <- function(df, response, group, group_levels, scenario_label) {
  pw <- pairwise_quiet(df, response, group)
  if (is.null(pw)) return(NULL)

  letters_tbl <- cld_from_pairwise(pw)
  if (is.null(letters_tbl)) return(NULL)

  # Headroom offset based on the overall data range so letters clear the
  # violin (which extends past the data max because trim = FALSE).
  rng <- range(df[[response]], na.rm = TRUE)
  offset <- diff(rng) * 0.18

  # Place every letter at one common height (overall max + offset) so they
  # form a level row above all violins. Per-group max positioning makes
  # low-n groups (e.g. category 2) sit lower and overlap their violin tail.
  y_common <- rng[2] + offset
  
  letters_tbl %>%
    mutate(y = y_common,
           scenario = scenario_label,
           group = factor(group, levels = group_levels))
}

# Violin plot (RCP 4.5 + 8.5 faceted) with CLD letters on top.
make_violin_cld <- function(df45, df85, group, group_levels,
                            title, subtitle, x_lab, file_name,
                            response = "trend",
                            dir = plot_dir) {
  long45 <- df45 %>%
    filter(!is.na(.data[[group]]), !is.na(.data[[response]])) %>%
    transmute(grp = factor(.data[[group]], levels = group_levels),
              scenario = "RCP 4.5", value = .data[[response]])
  long85 <- df85 %>%
    filter(!is.na(.data[[group]]), !is.na(.data[[response]])) %>%
    transmute(grp = factor(.data[[group]], levels = group_levels),
              scenario = "RCP 8.5", value = .data[[response]])
  pdata <- bind_rows(long45, long85)
  if (nrow(pdata) == 0) return(invisible(NULL))

  count_data <- pdata %>%
    group_by(grp, scenario) %>%
    summarise(n = n(), .groups = "drop")

  # CLD letters per scenario
  cld45 <- compute_cld_layer(df45, response, group, group_levels, "RCP 4.5")
  cld85 <- compute_cld_layer(df85, response, group, group_levels, "RCP 8.5")
  cld_data <- bind_rows(cld45, cld85)
  if (!is.null(cld_data) && nrow(cld_data) > 0) {
    cld_data <- cld_data %>% rename(grp = group)
  }

  p <- ggplot(pdata, aes(x = grp, y = value, fill = scenario)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white", alpha = 0.6) +
    geom_text(data = count_data, aes(x = grp, y = -Inf, label = n),
              vjust = -0.5, size = 2.5, inherit.aes = FALSE) +
    facet_wrap(~ scenario, ncol = 1) +
    scale_fill_manual(values = c("RCP 4.5" = "#0072B2",   # Okabe-Ito blue
                                 "RCP 8.5" = "#D55E00")) + # Okabe-Ito vermillion
    labs(title = title, subtitle = subtitle,
         x = x_lab, y = "Annual population trend (% per year)") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title    = element_text(size = 18, face = "bold"),
          plot.subtitle = element_text(size = 13),
          strip.text    = element_text(size = 14, face = "bold"),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          axis.text.x  = element_text(size = 12, face = "bold"),
          axis.text.y  = element_text(size = 12, face = "bold"))

  if (!is.null(cld_data) && nrow(cld_data) > 0) {
    p <- p + geom_text(data = cld_data,
                       aes(x = grp, y = y, label = cld),
                       vjust = 0, fontface = "bold", size = 4,
                       inherit.aes = FALSE)
  }

  out_f <- file.path(dir, file_name)
  ggsave(out_f, p, width = 8, height = 8, dpi = 150)
  cat("Saved plot:", basename(out_f), "\n")
  invisible(out_f)
}


# ==========================================================================
# PART 1: All 8 SDM categories (0-7) — all species ####
# ==========================================================================

stats_file_8cat <- file.path(stats_dir, paste0(land_cover, "_category_stats_8categories.txt"))
sink(stats_file_8cat)

cat("\n##################################################\n")
cat("# PART 1: All 8 SDM categories (0-7) — all species\n")
cat("#   Response: trend (Trend)\n")
cat("##################################################\n")

run_category_tests(analysis_45, "trend", "category", "RCP 4.5", "Trend")
run_category_tests(analysis_85, "trend", "category", "RCP 8.5", "Trend")

sink()
cat("Saved Part 1 (all-species 8-category) statistics to:", stats_file_8cat, "\n")

make_violin_cld(
  analysis_45, analysis_85,
  group = "category", group_levels = cat8_levels,
  title = paste0("All ", land_cover, " species combined"),
  subtitle = paste0(nrow(target_sdm), " route-species observations across ",
                    length(unique(target_sdm$route)), " unique routes and ",
                    length(unique(target_sdm$species_code)), " species"),
  x_lab = "SDM classified change category",
  file_name = paste0(land_cover, "_ALL_species_trend_by_rcp_cld.png")
)


# ==========================================================================
# PART 2: All 8 SDM categories (0-7) — per species ####
# ==========================================================================

stats_file_species_8cat <- file.path(stats_dir,
                                      paste0(land_cover, "_category_stats_by_species_8categories.txt"))
sink(stats_file_species_8cat)

cat("\n##################################################\n")
cat("# PART 2: All 8 SDM categories (0-7) — per species\n")
cat("#   Response: trend (Trend)\n")
cat("##################################################\n")

for (sp in target_species_list) {
  sp_name <- unique(target_sdm$species[target_sdm$species_code == sp])

  cat("\n==========================================================\n")
  cat(" ", sp_name, " (", sp, ")\n")
  cat("==========================================================\n")

  sp_45 <- analysis_45 %>% filter(species_code == sp)
  sp_85 <- analysis_85 %>% filter(species_code == sp)

  cat("\n--- RCP 4.5: All 8 categories ---\n")
  run_category_tests(sp_45, "trend", "category",
                     paste0(sp, " | RCP 4.5"), "Trend")

  cat("\n--- RCP 8.5: All 8 categories ---\n")
  run_category_tests(sp_85, "trend", "category",
                     paste0(sp, " | RCP 8.5"), "Trend")
}

sink()
cat("Saved Part 2 (per-species 8-category) statistics to:",
    stats_file_species_8cat, "\n")

for (sp in target_species_list) {
  sp_name <- unique(target_sdm$species[target_sdm$species_code == sp])
  sp_45 <- analysis_45 %>% filter(species_code == sp)
  sp_85 <- analysis_85 %>% filter(species_code == sp)

  make_violin_cld(
    sp_45, sp_85,
    group = "category", group_levels = cat8_levels,
    title = paste0(sp_name, " (", sp, ") (group:", land_cover, ")"),
    subtitle = paste0("n = ", nrow(sp_45), " (RCP4.5), ",
                      nrow(sp_85), " (RCP8.5) route observations"),
    x_lab = "SDM classified change category",
    file_name = paste0(land_cover, "_", sp, "_trend_by_rcp_cld.png"),
    dir = per_species_dir
  )
}


# ==========================================================================
# PART 3: Grouped categories (Contraction/Stable/Expansion) — all species ####
# ==========================================================================

stats_file_grp <- file.path(stats_dir, paste0(land_cover, "_category_stats_grouped.txt"))
sink(stats_file_grp)

cat("\n##################################################\n")
cat("# PART 3: Grouped categories — all species\n")
cat("#   Contraction = 1,2,3 | Stable = 4 | Expansion = 5,6,7\n")
cat("#   (category 0 = never suitable, excluded)\n")
cat("#   Response: trend (Trend)\n")
cat("##################################################\n")

run_category_tests(analysis_45_grp, "trend", "change_group", "RCP 4.5", "Trend")
run_category_tests(analysis_85_grp, "trend", "change_group", "RCP 8.5", "Trend")

sink()
cat("Saved Part 3 (all-species grouped) statistics to:", stats_file_grp, "\n")

make_violin_cld(
  analysis_45_grp, analysis_85_grp,
  group = "change_group", group_levels = grp_levels,
  title = paste0("All ", land_cover, " species combined"),
  subtitle = paste0(nrow(target_sdm), " route-species observations across ",
                    length(unique(target_sdm$route)), " unique routes and ",
                    length(unique(target_sdm$species_code)), " species"),
  x_lab = "Range change group",
  file_name = paste0(land_cover, "_ALL_species_trend_by_group_cld.png")
)


# ==========================================================================
# PART 4: Grouped categories (Contraction/Stable/Expansion) — per species ####
# ==========================================================================

stats_file_species_grp <- file.path(stats_dir,
                                     paste0(land_cover, "_category_stats_by_species_grouped.txt"))
sink(stats_file_species_grp)

cat("\n##################################################\n")
cat("# PART 4: Grouped categories — per species\n")
cat("#   Contraction = 1,2,3 | Stable = 4 | Expansion = 5,6,7\n")
cat("#   Response: trend (Trend)\n")
cat("##################################################\n")

for (sp in target_species_list) {
  sp_name <- unique(target_sdm$species[target_sdm$species_code == sp])

  cat("\n==========================================================\n")
  cat(" ", sp_name, " (", sp, ")\n")
  cat("==========================================================\n")

  sp_45 <- analysis_45_grp %>% filter(species_code == sp)
  sp_85 <- analysis_85_grp %>% filter(species_code == sp)

  cat("\n--- RCP 4.5: Grouped categories ---\n")
  run_category_tests(sp_45, "trend", "change_group",
                     paste0(sp, " | RCP 4.5"), "Trend")

  cat("\n--- RCP 8.5: Grouped categories ---\n")
  run_category_tests(sp_85, "trend", "change_group",
                     paste0(sp, " | RCP 8.5"), "Trend")
}

sink()
cat("Saved Part 4 (per-species grouped) statistics to:",
    stats_file_species_grp, "\n")

for (sp in target_species_list) {
  sp_name <- unique(target_sdm$species[target_sdm$species_code == sp])
  sp_45 <- analysis_45_grp %>% filter(species_code == sp)
  sp_85 <- analysis_85_grp %>% filter(species_code == sp)

  make_violin_cld(
    sp_45, sp_85,
    group = "change_group", group_levels = grp_levels,
    title = paste0(sp_name, " (", sp, ") (group:", land_cover, ")"),
    subtitle = paste0("n = ", nrow(sp_45), " (RCP4.5), ",
                      nrow(sp_85), " (RCP8.5) route observations"),
    x_lab = "Range change group",
    file_name = paste0(land_cover, "_", sp, "_trend_by_group_cld.png"),
    dir = per_species_dir
  )
}


# ==========================================================================
# PART 5: ALL land covers combined (ignore land_cover filter) ####
#   Pools every *_route_trends_sdm.csv in out_dir, across every land-cover
#   group present in spp_names_codes_group_aou.csv (e.g. grasslands and
#   aridlands), into one dataset for the all-species 8-category and grouped
#   analyses + plots.
# ==========================================================================

cat("\nCombining all species SDM files across ALL land covers...\n")

# Reuse the already-loaded, already-joined data (all_sdm_raw carries every
# species' Group from the lookup table) instead of re-reading the files.
if (nrow(all_sdm_raw) == 0) {
  warning("No SDM CSVs found in ", out_dir, " for the all-land-cover analysis.")
} else {
  all_sdm <- all_sdm_raw %>% rename(land_cover = Group)

  land_covers_used <- sort(unique(all_sdm$land_cover))
  cat("Land covers:", paste(land_covers_used, collapse = ", "), "\n")
  cat("Total rows:", nrow(all_sdm), "\n")
  cat("Species:", length(unique(all_sdm$species_code)), "\n")

  # Build combined analysis data (mirrors the per-land-cover setup above).
  analysis_45_all <- all_sdm %>%
    filter(!is.na(rcp45)) %>%
    transmute(route, species_code, category = rcp45, trend,
              route_num = as.integer(sub("^\\d+-", "", route)))

  analysis_85_all <- all_sdm %>%
    filter(!is.na(rcp85)) %>%
    transmute(route, species_code, category = rcp85, trend,
              route_num = as.integer(sub("^\\d+-", "", route)))

  analysis_45_all_grp <- analysis_45_all %>%
    mutate(change_group = factor(group_category(category),
                                 levels = c("Contraction", "Stable", "Expansion")))

  analysis_85_all_grp <- analysis_85_all %>%
    mutate(change_group = factor(group_category(category),
                                 levels = c("Contraction", "Stable", "Expansion")))

  lc_tag <- "all_lc"

  # ----- 8 SDM categories (0-7) — all species, all land covers -----
  stats_file_all_8cat <- file.path(stats_dir,
                                    paste0(lc_tag, "_category_stats_8categories.txt"))
  sink(stats_file_all_8cat)

  cat("\n##################################################\n")
  cat("# PART 5: All 8 SDM categories (0-7) — all species, ALL land covers\n")
  cat("#   Land covers pooled:", paste(land_covers_used, collapse = ", "), "\n")
  cat("#   Response: trend (Trend)\n")
  cat("##################################################\n")

  run_category_tests(analysis_45_all, "trend", "category", "RCP 4.5", "Trend")
  run_category_tests(analysis_85_all, "trend", "category", "RCP 8.5", "Trend")

  sink()
  cat("Saved Part 5 (all-land-cover 8-category) statistics to:", stats_file_all_8cat, "\n")

  make_violin_cld(
    analysis_45_all, analysis_85_all,
    group = "category", group_levels = cat8_levels,
    title = "All species combined — all land covers",
    subtitle = paste0(nrow(all_sdm), " route-species observations across ",
                      length(unique(all_sdm$route)), " unique routes, ",
                      length(unique(all_sdm$species_code)), " species, and ",
                      length(land_covers_used), " land covers"),
    x_lab = "SDM classified change category",
    file_name = paste0(lc_tag, "_ALL_species_trend_by_rcp_cld.png")
  )

  # ----- Grouped categories — all species, all land covers -----
  stats_file_all_grp <- file.path(stats_dir,
                                   paste0(lc_tag, "_category_stats_grouped.txt"))
  sink(stats_file_all_grp)

  cat("\n##################################################\n")
  cat("# PART 5: Grouped categories — all species, ALL land covers\n")
  cat("#   Contraction = 1,2,3 | Stable = 4 | Expansion = 5,6,7\n")
  cat("#   (category 0 = never suitable, excluded)\n")
  cat("#   Land covers pooled:", paste(land_covers_used, collapse = ", "), "\n")
  cat("#   Response: trend (Trend)\n")
  cat("##################################################\n")

  run_category_tests(analysis_45_all_grp, "trend", "change_group", "RCP 4.5", "Trend")
  run_category_tests(analysis_85_all_grp, "trend", "change_group", "RCP 8.5", "Trend")

  sink()
  cat("Saved Part 5 (all-land-cover grouped) statistics to:", stats_file_all_grp, "\n")

  make_violin_cld(
    analysis_45_all_grp, analysis_85_all_grp,
    group = "change_group", group_levels = grp_levels,
    title = "All species combined — all land covers",
    subtitle = paste0(nrow(all_sdm), " route-species observations across ",
                      length(unique(all_sdm$route)), " unique routes, ",
                      length(unique(all_sdm$species_code)), " species, and ",
                      length(land_covers_used), " land covers"),
    x_lab = "Range change group",
    file_name = paste0(lc_tag, "_ALL_species_trend_by_group_cld.png")
  )
}

# ==========================================================================
# PART 6: Tier A species only — combined visualization ####
#   Subsets the target land_cover dataset to the 13 Tier-A (converged,
#   reliable) species and reruns the all-species-combined stats + plots for
#   both the 8-category and grouped views.
# ==========================================================================

# Tier A species: the 13 grassland species that converged cleanly
# (R-hat <= 1.02, ESS at/near target) per
# output/model_diagnostics/grassland_diagnostics_SUMMARY.md. The combined
# Tier-A visualizations (Part 6) pool only these reliable fits.
tierA_species <- c("BOSP", "LEPC", "WTHA", "NESP", "SPPI", "MCLO", "MOPL",
                   "SEOW", "STGR", "LCSP", "GRPC", "HESP", "GRAP")

cat("\nBuilding Tier-A-only combined analysis...\n")

tierA_present <- intersect(tierA_species, unique(target_sdm$species_code))
tierA_missing <- setdiff(tierA_species, unique(target_sdm$species_code))
if (length(tierA_missing) > 0) {
  cat("Tier-A species with no SDM CSV (skipped):",
      paste(tierA_missing, collapse = ", "), "\n")
}

if (length(tierA_present) == 0) {
  warning("No Tier-A species found in ", land_cover,
          " SDM data; skipping Part 6.")
} else {
  tierA_sdm <- target_sdm %>% filter(species_code %in% tierA_present)

  cat("Tier-A species included:", length(tierA_present), "\n")
  cat("Tier-A total rows:", nrow(tierA_sdm), "\n")

  analysis_45_A     <- analysis_45     %>% filter(species_code %in% tierA_present)
  analysis_85_A     <- analysis_85     %>% filter(species_code %in% tierA_present)
  analysis_45_A_grp <- analysis_45_grp %>% filter(species_code %in% tierA_present)
  analysis_85_A_grp <- analysis_85_grp %>% filter(species_code %in% tierA_present)

  tierA_tag <- paste0(land_cover, "_tierA")

  # ----- 8 SDM categories (0-7) — Tier A species -----
  stats_file_A_8cat <- file.path(stats_dir,
                                 paste0(tierA_tag, "_category_stats_8categories.txt"))
  sink(stats_file_A_8cat)

  cat("\n##################################################\n")
  cat("# PART 6: All 8 SDM categories (0-7) — Tier A species only\n")
  cat("#   Tier-A species:", paste(tierA_present, collapse = ", "), "\n")
  cat("#   Response: trend (Trend)\n")
  cat("##################################################\n")

  run_category_tests(analysis_45_A, "trend", "category", "RCP 4.5", "Trend")
  run_category_tests(analysis_85_A, "trend", "category", "RCP 8.5", "Trend")

  sink()
  cat("Saved Part 6 (Tier-A 8-category) statistics to:", stats_file_A_8cat, "\n")

  make_violin_cld(
    analysis_45_A, analysis_85_A,
    group = "category", group_levels = cat8_levels,
    title = paste0("Tier A ", land_cover, " species combined"),
    subtitle = paste0(nrow(tierA_sdm), " route-species observations across ",
                      length(unique(tierA_sdm$route)), " unique routes and ",
                      length(tierA_present), " converged species"),
    x_lab = "SDM classified change category",
    file_name = paste0(tierA_tag, "_ALL_species_trend_by_rcp_cld.png")
  )

  # ----- Grouped categories — Tier A species -----
  stats_file_A_grp <- file.path(stats_dir,
                                paste0(tierA_tag, "_category_stats_grouped.txt"))
  sink(stats_file_A_grp)

  cat("\n##################################################\n")
  cat("# PART 6: Grouped categories — Tier A species only\n")
  cat("#   Contraction = 1,2,3 | Stable = 4 | Expansion = 5,6,7\n")
  cat("#   (category 0 = never suitable, excluded)\n")
  cat("#   Tier-A species:", paste(tierA_present, collapse = ", "), "\n")
  cat("#   Response: trend (Trend)\n")
  cat("##################################################\n")

  run_category_tests(analysis_45_A_grp, "trend", "change_group", "RCP 4.5", "Trend")
  run_category_tests(analysis_85_A_grp, "trend", "change_group", "RCP 8.5", "Trend")

  sink()
  cat("Saved Part 6 (Tier-A grouped) statistics to:", stats_file_A_grp, "\n")

  make_violin_cld(
    analysis_45_A_grp, analysis_85_A_grp,
    group = "change_group", group_levels = grp_levels,
    title = paste0("Tier A ", land_cover, " species combined"),
    subtitle = paste0(nrow(tierA_sdm), " route-species observations across ",
                      length(unique(tierA_sdm$route)), " unique routes and ",
                      length(tierA_present), " converged species"),
    x_lab = "Range change group",
    file_name = paste0(tierA_tag, "_ALL_species_trend_by_group_cld.png")
  )
}


cat("\n=== Statistical analysis + visualization complete ===\n")
