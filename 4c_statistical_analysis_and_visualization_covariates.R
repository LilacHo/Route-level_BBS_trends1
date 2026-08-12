# 4c_statistical_analysis_and_visualization_covariates.R
#
# Covariate-pipeline counterpart of 4_statistical_analysis_and_visualization.R
# — identical logic, different input/output directories, so it can run
# against 3c_add_SDM_covariates.R's output without touching
# 4_statistical_analysis_and_visualization.R or colliding with its
# output/species_routes_sdm_stats / output/species_routes_sdm_plot
# directories (which belong to the original, non-covariate pipeline).
#
# Combines the route-level trend statistical analysis with violin plots
# annotated with compact letter displays (CLD).
#
# The CLD letters are derived from the same pairwise Wilcoxon (BH-adjusted)
# tests used in the statistical analysis: groups that share a letter are NOT
# significantly different (alpha = 0.05).
#
# Focus response throughout: trend (annual % change in route-level counts,
# written by 2c_generate_route_trend_csvs_covariates.R's per-species-per-model
# CSVs and carried through 3c_add_SDM_covariates.R).
#
# Every input file has a "model" column (base / anthro), since 3c's input
# files are per-species-PER-MODEL. This script analyzes exactly ONE
# model_tag per run (set below) — all rows for every other model are
# filtered out immediately after reading, before any test or plot is built,
# so nothing gets pooled across models. To compare models, rerun this
# script once per model_tag (edit the setting, or set `model_tag <- "..."`
# before sourcing this file) — every stats/plot output filename below
# includes model_tag, so separate runs don't overwrite each other's output.
#
# Since 1c_species_iCAR_covariates.R / 2c now fit and combine EVERY species
# across all 12 groups in one run (not one bird_group at a time), each
# input file also carries a "group" column (written directly by 2c). Set
# bird_group below to one real group (e.g. "aridlands") to reproduce the
# original per-group Parts 1-4 breakdown, or set it to NA to skip the
# group filter and pool every species/group together (PART 5 below, which
# used to be an optional add-on, is now effectively the default view given
# 1c's all-species scope, and runs automatically when bird_group is NA).

library(here)
library(tidyverse)
library(multcompView) # multcompView converts pairwise p-values into compact letter displays.

here::i_am("4c_statistical_analysis_and_visualization_covariates.R")

# Settings -----------------------------------------------------------------
bird_group <- NA   # one of the 12 Group values in
                   # data/spp_names_codes_group_aou.csv (e.g. "aridlands",
                   # "boreal_forests", ...) to restrict Parts 1-4 to that
                   # group, or NA to skip Parts 1-4 entirely and only run
                   # Part 5 (all groups pooled) -- see note above.

model_tag <- "base" # Change as needed ("base" or "anthro")

out_dir   <- here::here("output", "species_routes_covariates", "per_species_sdm")
stats_dir <- here::here("output", "species_routes_covariates", "per_species_sdm_stats")
plot_dir  <- here::here("output", "species_routes_covariates", "per_species_sdm_plot")
per_species_dir <- file.path(plot_dir, "per_species")
if (!dir.exists(stats_dir))       dir.create(stats_dir,       recursive = TRUE)
if (!dir.exists(plot_dir))        dir.create(plot_dir,        recursive = TRUE)
if (!dir.exists(per_species_dir)) dir.create(per_species_dir, recursive = TRUE)

# Read every SDM CSV written by 3c_add_SDM_covariates.R. Each row already
# carries its own "group" column (written directly by
# 2c_generate_route_trend_csvs_covariates.R and passed through unchanged by
# 3c), so no separate species lookup/join is needed here.
sdm_files <- list.files(out_dir, pattern = "_route_trends_sdm\\.csv$",
                        full.names = TRUE)

if (length(sdm_files) == 0) {
  stop("No SDM CSVs found in ", out_dir, " — run 3c_add_SDM_covariates.R first.")
}

cat("\nCombining all species SDM files...\n")

all_sdm_raw_unfiltered <- sdm_files %>% map_dfr(read.csv)

if (!"model" %in% names(all_sdm_raw_unfiltered)) {
  stop("No 'model' column found in the SDM CSVs in ", out_dir, " — expected ",
       "one per-species-per-model file per 2c_generate_route_trend_csvs_covariates.R's ",
       "output. Re-run 2c/3c if these files predate that per-model split.")
}
if (!"group" %in% names(all_sdm_raw_unfiltered)) {
  stop("No 'group' column found in the SDM CSVs in ", out_dir, " — these files ",
       "predate 2c writing 'group' directly; re-run 2c_generate_route_trend_csvs_covariates.R ",
       "and 3c_add_SDM_covariates.R to regenerate them.")
}

models_available <- sort(unique(all_sdm_raw_unfiltered$model))
if (!model_tag %in% models_available) {
  stop("model_tag = '", model_tag, "' not found in ", out_dir,
       " — models available: ", paste(models_available, collapse = ", "),
       ". Set model_tag to one of those (or run 2c/3c for the model you want first).")
}

# Restrict to exactly one model BEFORE any downstream filtering/analysis, so
# every stats/plot output below reflects this one model only — nothing gets
# pooled across models.
all_sdm_raw <- all_sdm_raw_unfiltered %>% filter(model == model_tag)

cat("model_tag:", model_tag, "(", nrow(all_sdm_raw), "of", nrow(all_sdm_raw_unfiltered),
    "rows across all models present in", out_dir, ")\n")

# bird_group == NA means "pool every group" (see settings comment above);
# otherwise restrict to that one group, same as before. run_label is used
# everywhere below that used to use bird_group directly for file/title
# text, so filenames never end up with a literal "NA" in them.
if (is.na(bird_group)) {
  target_sdm <- all_sdm_raw
  run_label  <- "all_groups"
  cat("bird_group = NA -> pooling all groups\n")
} else {
  target_sdm <- all_sdm_raw %>% filter(group == bird_group)
  run_label  <- bird_group
}

if (nrow(target_sdm) == 0) {
  stop("No SDM rows matched bird_group = '", bird_group, "' and model_tag = '",
       model_tag, "'. Groups present: ", paste(sort(unique(all_sdm_raw$group)), collapse = ", "))
}

cat("Total rows:", nrow(target_sdm), "\n")
cat(run_label, "species:", length(unique(target_sdm$species_code)), "\n")
if (is.na(bird_group)) {
  cat("Groups pooled:", paste(sort(unique(target_sdm$group)), collapse = ", "), "\n")
}

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

stats_file_8cat <- file.path(stats_dir, paste0(run_label, "_", model_tag, "_category_stats_8categories.txt"))
sink(stats_file_8cat)

cat("\n##################################################\n")
cat("# PART 1: All 8 SDM categories (0-7) — all species\n")
cat("#   Model:", model_tag, "\n")
cat("#   Response: trend (Trend)\n")
cat("##################################################\n")

run_category_tests(analysis_45, "trend", "category", "RCP 4.5", "Trend")
run_category_tests(analysis_85, "trend", "category", "RCP 8.5", "Trend")

sink()
cat("Saved Part 1 (all-species 8-category) statistics to:", stats_file_8cat, "\n")

make_violin_cld(
  analysis_45, analysis_85,
  group = "category", group_levels = cat8_levels,
  title = paste0("All ", run_label, " species combined (", model_tag, ")"),
  subtitle = paste0(nrow(target_sdm), " route-species observations across ",
                    length(unique(target_sdm$route)), " unique routes and ",
                    length(unique(target_sdm$species_code)), " species"),
  x_lab = "SDM classified change category",
  file_name = paste0(run_label, "_", model_tag, "_ALL_species_trend_by_rcp_cld.png")
)


# ==========================================================================
# PART 2: All 8 SDM categories (0-7) — per species ####
# ==========================================================================

stats_file_species_8cat <- file.path(stats_dir,
                                      paste0(run_label, "_", model_tag, "_category_stats_by_species_8categories.txt"))
sink(stats_file_species_8cat)

cat("\n##################################################\n")
cat("# PART 2: All 8 SDM categories (0-7) — per species\n")
cat("#   Model:", model_tag, "\n")
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
  sp_grp  <- unique(target_sdm$group[target_sdm$species_code == sp])[1]
  sp_45 <- analysis_45 %>% filter(species_code == sp)
  sp_85 <- analysis_85 %>% filter(species_code == sp)

  make_violin_cld(
    sp_45, sp_85,
    group = "category", group_levels = cat8_levels,
    title = paste0(sp_name, " (", sp, ") (group:", sp_grp, ", model:", model_tag, ")"),
    subtitle = paste0("n = ", nrow(sp_45), " (RCP4.5), ",
                      nrow(sp_85), " (RCP8.5) route observations"),
    x_lab = "SDM classified change category",
    file_name = paste0(run_label, "_", model_tag, "_", sp, "_trend_by_rcp_cld.png"),
    dir = per_species_dir
  )
}


# ==========================================================================
# PART 3: Grouped categories (Contraction/Stable/Expansion) — all species ####
# ==========================================================================

stats_file_grp <- file.path(stats_dir, paste0(run_label, "_", model_tag, "_category_stats_grouped.txt"))
sink(stats_file_grp)

cat("\n##################################################\n")
cat("# PART 3: Grouped categories — all species\n")
cat("#   Contraction = 1,2,3 | Stable = 4 | Expansion = 5,6,7\n")
cat("#   (category 0 = never suitable, excluded)\n")
cat("#   Model:", model_tag, "\n")
cat("#   Response: trend (Trend)\n")
cat("##################################################\n")

run_category_tests(analysis_45_grp, "trend", "change_group", "RCP 4.5", "Trend")
run_category_tests(analysis_85_grp, "trend", "change_group", "RCP 8.5", "Trend")

sink()
cat("Saved Part 3 (all-species grouped) statistics to:", stats_file_grp, "\n")

make_violin_cld(
  analysis_45_grp, analysis_85_grp,
  group = "change_group", group_levels = grp_levels,
  title = paste0("All ", run_label, " species combined (", model_tag, ")"),
  subtitle = paste0(nrow(target_sdm), " route-species observations across ",
                    length(unique(target_sdm$route)), " unique routes and ",
                    length(unique(target_sdm$species_code)), " species"),
  x_lab = "Range change group",
  file_name = paste0(run_label, "_", model_tag, "_ALL_species_trend_by_group_cld.png")
)


# ==========================================================================
# PART 4: Grouped categories (Contraction/Stable/Expansion) — per species ####
# ==========================================================================

stats_file_species_grp <- file.path(stats_dir,
                                     paste0(run_label, "_", model_tag, "_category_stats_by_species_grouped.txt"))
sink(stats_file_species_grp)

cat("\n##################################################\n")
cat("# PART 4: Grouped categories — per species\n")
cat("#   Contraction = 1,2,3 | Stable = 4 | Expansion = 5,6,7\n")
cat("#   Model:", model_tag, "\n")
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
  sp_grp  <- unique(target_sdm$group[target_sdm$species_code == sp])[1]
  sp_45 <- analysis_45_grp %>% filter(species_code == sp)
  sp_85 <- analysis_85_grp %>% filter(species_code == sp)

  make_violin_cld(
    sp_45, sp_85,
    group = "change_group", group_levels = grp_levels,
    title = paste0(sp_name, " (", sp, ") (group:", sp_grp, ", model:", model_tag, ")"),
    subtitle = paste0("n = ", nrow(sp_45), " (RCP4.5), ",
                      nrow(sp_85), " (RCP8.5) route observations"),
    x_lab = "Range change group",
    file_name = paste0(run_label, "_", model_tag, "_", sp, "_trend_by_group_cld.png"),
    dir = per_species_dir
  )
}


# PART 5 (all land covers combined) used to live here as a separate,
# commented-out block. It's no longer needed as a separate section: Parts
# 1-4 above now do the same thing natively whenever bird_group <- NA (see
# the settings block and run_label/target_sdm construction near the top of
# this script) — target_sdm already becomes "every group pooled" and every
# stats/plot filename already uses run_label ("all_groups" in that case), so
# there's nothing left for a separate Part 5 to add. Set bird_group to a
# specific group instead of NA to get the original single-group behavior.

cat("\n=== Statistical analysis + visualization complete (covariate pipeline, model_tag = ",
    model_tag, ", bird_group = ", run_label, ") ===\n", sep = "")
