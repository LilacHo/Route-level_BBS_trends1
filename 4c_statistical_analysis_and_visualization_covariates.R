# 4c_statistical_analysis_and_visualization_covariates.R
#
# Statistical comparison and violin plots of route-level trend across SDM
# classified-change categories, with compact letter displays (CLD) from
# pairwise Wilcoxon tests (BH-adjusted; groups sharing a letter are not
# significantly different at alpha = 0.05).
#
# Analyzes exactly ONE model_tag per run (set below) -- rerun once per tag
# to compare base / grassland_habitat / grassland_anthro /
# grassland_habitat_to_anthro; every output filename includes model_tag so
# separate runs don't overwrite each other.
#
# Input: output/species_routes_covariates/per_species_sdm/*.csv (from 3c),
# data/spp_names_codes_group_aou.csv.
#
# Output: stats reports (.txt) in output/species_routes_covariates/per_species_sdm_stats/,
# violin plots (.png, all-species and per-species) in
# output/species_routes_covariates/per_species_sdm_plot/, covering both the
# raw 8 SDM categories (0-7) and a grouped Contraction/Stable/Expansion view.

library(here)
library(tidyverse)
library(multcompView) # converts pairwise p-values into compact letter displays

here::i_am("4c_statistical_analysis_and_visualization_covariates.R")

# Settings -----------------------------------------------------------------
land_cover <- "grasslands"   # species are matched to this group via
                             # data/spp_names_codes_group_aou.csv

model_tag <- "base" # Change as needed
# (base / grassland_habitat / grassland_anthro / grassland_habitat_to_anthro)

out_dir   <- here::here("output", "species_routes_covariates", "per_species_sdm")
stats_dir <- here::here("output", "species_routes_covariates", "per_species_sdm_stats")
plot_dir  <- here::here("output", "species_routes_covariates", "per_species_sdm_plot")
per_species_dir <- file.path(plot_dir, "per_species")
if (!dir.exists(stats_dir))       dir.create(stats_dir,       recursive = TRUE)
if (!dir.exists(plot_dir))        dir.create(plot_dir,        recursive = TRUE)
if (!dir.exists(per_species_dir)) dir.create(per_species_dir, recursive = TRUE)

# Species lookup (Code -> Group), to filter the SDM CSVs by land_cover ------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv")) %>%
  distinct(Code, Group)

sdm_files <- list.files(out_dir, pattern = "_route_trends_sdm\\.csv$",
                        full.names = TRUE)

if (length(sdm_files) == 0) {
  stop("No SDM CSVs found in ", out_dir, " — run 3c_add_SDM_covariates.R first.")
}

cat("\nCombining all species SDM files...\n")

all_sdm_raw_unfiltered <- sdm_files %>%
  map_dfr(read.csv) %>%
  left_join(spp_df, by = c("species_code" = "Code"))

if (!"model" %in% names(all_sdm_raw_unfiltered)) {
  stop("No 'model' column found in the SDM CSVs in ", out_dir, " — expected ",
       "one per-species-per-model file per 2c_generate_route_trend_csvs_covariates.R's ",
       "output. Re-run 2c/3c if these files predate that per-model split.")
}

models_available <- sort(unique(all_sdm_raw_unfiltered$model))
if (!model_tag %in% models_available) {
  stop("model_tag = '", model_tag, "' not found in ", out_dir,
       " — models available: ", paste(models_available, collapse = ", "),
       ". Set model_tag to one of those (or run 2c/3c for the model you want first).")
}

# Restrict to one model before any filtering/analysis, so nothing pools
# across models.
all_sdm_raw <- all_sdm_raw_unfiltered %>% filter(model == model_tag)

cat("model_tag:", model_tag, "(", nrow(all_sdm_raw), "of", nrow(all_sdm_raw_unfiltered),
    "rows across all models present in", out_dir, ")\n")

target_sdm <- all_sdm_raw %>% filter(Group == land_cover)

if (nrow(target_sdm) == 0) {
  stop("No SDM rows matched land_cover = '", land_cover, "' and model_tag = '",
       model_tag, "' after joining species groups from spp_names_codes_group_aou.csv")
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

# Runs Kruskal-Wallis + pairwise Wilcoxon (BH-corrected), prints a report,
# and returns the fitted objects so the pairwise matrix can feed the CLD.
run_category_tests <- function(df, response, group, scenario_label, response_label) {
  d <- df %>%
    filter(!is.na(.data[[response]]), !is.na(.data[[group]])) %>%
    mutate(grp = factor(.data[[group]]))

  cat("\n--------------------------------------------------\n")
  cat(scenario_label, "|", response_label, "by", group, "\n")
  cat("--------------------------------------------------\n")

  if (length(unique(d$grp)) < 2) {
    cat("  Not enough groups for testing.\n")
    return(invisible(NULL))
  }

  summ_tab <- d %>%
    group_by(grp) %>%
    summarise(n = n(),
              median = median(.data[[response]]),
              mean = mean(.data[[response]]),
              .groups = "drop")
  cat("\nGroup summary:\n")
  print(as.data.frame(summ_tab), row.names = FALSE)

  kw <- kruskal.test(d[[response]] ~ d$grp)
  cat("\nKruskal-Wallis: chi-sq =", round(kw$statistic, 2),
      ", df =", kw$parameter,
      ", p =", format.pval(kw$p.value, digits = 3), "\n")

  pw <- pairwise.wilcox.test(d[[response]], d$grp, p.adjust.method = "BH")
  cat("\nPairwise Wilcoxon (BH-adjusted p-values):\n")
  print(round(pw$p.value, 4))

  invisible(list(kruskal = kw, pairwise = pw, data = d))
}

# Builds a full symmetric p-value matrix from a pairwise.wilcox.test result
# (which returns lower-triangular) since multcompLetters needs a named
# vector of p-values for all pairs.
cld_from_pairwise <- function(pw, alpha = 0.05) {
  if (is.null(pw) || is.null(pw$p.value)) return(NULL)

  pmat <- pw$p.value
  groups <- union(rownames(pmat), colnames(pmat))

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

# Quiet pairwise Wilcoxon (BH) for CLD purposes only -- no printing, so the
# plotting layer doesn't duplicate the stats report already written to file.
pairwise_quiet <- function(df, response, group) {
  d <- df %>%
    filter(!is.na(.data[[response]]), !is.na(.data[[group]])) %>%
    mutate(grp = factor(.data[[group]]))
  if (length(unique(d$grp)) < 2) return(NULL)
  pairwise.wilcox.test(d[[response]], d$grp, p.adjust.method = "BH")
}

# Per-group CLD labels for one scenario/grouping: group, scenario, letter,
# and a shared y position above all violins (so low-n groups don't overlap
# their own violin tail).
compute_cld_layer <- function(df, response, group, group_levels, scenario_label) {
  pw <- pairwise_quiet(df, response, group)
  if (is.null(pw)) return(NULL)

  letters_tbl <- cld_from_pairwise(pw)
  if (is.null(letters_tbl)) return(NULL)

  rng <- range(df[[response]], na.rm = TRUE)
  offset <- diff(rng) * 0.18
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

stats_file_8cat <- file.path(stats_dir, paste0(land_cover, "_", model_tag, "_category_stats_8categories.txt"))
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
  title = paste0("All ", land_cover, " species combined (", model_tag, ")"),
  subtitle = paste0(nrow(target_sdm), " route-species observations across ",
                    length(unique(target_sdm$route)), " unique routes and ",
                    length(unique(target_sdm$species_code)), " species"),
  x_lab = "SDM classified change category",
  file_name = paste0(land_cover, "_", model_tag, "_ALL_species_trend_by_rcp_cld.png")
)


# ==========================================================================
# PART 2: All 8 SDM categories (0-7) — per species ####
# ==========================================================================

stats_file_species_8cat <- file.path(stats_dir,
                                      paste0(land_cover, "_", model_tag, "_category_stats_by_species_8categories.txt"))
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
  sp_45 <- analysis_45 %>% filter(species_code == sp)
  sp_85 <- analysis_85 %>% filter(species_code == sp)

  make_violin_cld(
    sp_45, sp_85,
    group = "category", group_levels = cat8_levels,
    title = paste0(sp_name, " (", sp, ") (group:", land_cover, ", model:", model_tag, ")"),
    subtitle = paste0("n = ", nrow(sp_45), " (RCP4.5), ",
                      nrow(sp_85), " (RCP8.5) route observations"),
    x_lab = "SDM classified change category",
    file_name = paste0(land_cover, "_", model_tag, "_", sp, "_trend_by_rcp_cld.png"),
    dir = per_species_dir
  )
}


# ==========================================================================
# PART 3: Grouped categories (Contraction/Stable/Expansion) — all species ####
# ==========================================================================

stats_file_grp <- file.path(stats_dir, paste0(land_cover, "_", model_tag, "_category_stats_grouped.txt"))
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
  title = paste0("All ", land_cover, " species combined (", model_tag, ")"),
  subtitle = paste0(nrow(target_sdm), " route-species observations across ",
                    length(unique(target_sdm$route)), " unique routes and ",
                    length(unique(target_sdm$species_code)), " species"),
  x_lab = "Range change group",
  file_name = paste0(land_cover, "_", model_tag, "_ALL_species_trend_by_group_cld.png")
)


# ==========================================================================
# PART 4: Grouped categories (Contraction/Stable/Expansion) — per species ####
# ==========================================================================

stats_file_species_grp <- file.path(stats_dir,
                                     paste0(land_cover, "_", model_tag, "_category_stats_by_species_grouped.txt"))
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
  sp_45 <- analysis_45_grp %>% filter(species_code == sp)
  sp_85 <- analysis_85_grp %>% filter(species_code == sp)

  make_violin_cld(
    sp_45, sp_85,
    group = "change_group", group_levels = grp_levels,
    title = paste0(sp_name, " (", sp, ") (group:", land_cover, ", model:", model_tag, ")"),
    subtitle = paste0("n = ", nrow(sp_45), " (RCP4.5), ",
                      nrow(sp_85), " (RCP8.5) route observations"),
    x_lab = "Range change group",
    file_name = paste0(land_cover, "_", model_tag, "_", sp, "_trend_by_group_cld.png"),
    dir = per_species_dir
  )
}


cat("\n=== Statistical analysis + visualization complete (covariate pipeline, model_tag = ",
    model_tag, ") ===\n", sep = "")
