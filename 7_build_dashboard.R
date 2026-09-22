# 7_build_dashboard.R
#
# Builds a static, self-contained local dashboard (dashboard/) for browsing,
# per species x model (base/anthro) x climate scenario (RCP4.5/RCP8.5):
#   - the by-RCP violin plot (8 SDM categories, both scenarios in one image)
#   - the by-group violin plot (Contraction/Stable/Expansion, both scenarios)
#   - the range-shift map (changes with model AND scenario)
#   - a Contraction/Stable/Expansion table (n, mean, sd) for both models
#     side by side
#
# Scope: only the 469 species that actually have model output (i.e. appear
# in output/species_routes_covariates/per_species_sdm/) -- this is the
# authoritative species list, taken straight from those CSVs rather than
# data/spp_names_codes_group_aou.csv (which has ~611 rows, most without any
# model output).
#
# No fetch()/server needed: manifest + table data are written as plain .js
# files (`const X = {...}`) loaded via <script src>, so dashboard/index.html
# works by double-clicking it, same as opening any other local file.
#
# Reads:
#   output/species_routes_covariates/per_species_sdm/*_base_route_trends_sdm.csv
#     (species list: code, name, group, slug -- one row read per file)
#   output/species_routes_covariates/per_species_sdm_plot/per_species/*.png
#     (by-RCP and by-group violins)
#   output/species_maps/<model>/<rcp>/*.png
#     (range-shift maps)
#   output/files/all_species_habitat_contraction_stable_expansion_2010_2025_{base,anthro}.xlsx
#     (Contraction/Stable/Expansion n/mean/sd table, both RCP sheets)
#
# Writes:
#   dashboard/data/manifest.js, dashboard/data/tables.js
#   dashboard/img/violin_rcp/<code>_<model>.webp
#   dashboard/img/violin_group/<code>_<model>.webp   (skipped where no source exists)
#   dashboard/img/maps/<code>_<model>_<rcp>.webp
#
# Idempotent: re-running overwrites all generated files; nothing here is
# resumable/skippable since a full rebuild only takes a few minutes.

library(here)
library(tidyverse)
library(openxlsx)
library(magick)
library(jsonlite)

here::i_am("7_build_dashboard.R")

# Settings -------------------------------------------------------------------
models    <- c("base", "anthro")
scenarios <- c("rcp45", "rcp85")
scenario_sheet <- c(rcp45 = "RCP4.5", rcp85 = "RCP8.5")
grp_levels <- c("Contraction", "Stable", "Expansion")

violin_width  <- 700   # px, violins are square (ggsave width=height=8in)
map_width     <- 1400  # px, maps are all fig_width=10in wide -> uniform scale factor
webp_quality  <- 82

sdm_dir      <- here::here("output", "species_routes_covariates", "per_species_sdm")
violin_dir   <- here::here("output", "species_routes_covariates", "per_species_sdm_plot", "per_species")
maps_src_dir <- here::here("output", "species_maps")
files_dir    <- here::here("output", "files")

dash_dir       <- here::here("dashboard")
data_dir       <- file.path(dash_dir, "data")
violin_rcp_out <- file.path(dash_dir, "img", "violin_rcp")
violin_grp_out <- file.path(dash_dir, "img", "violin_group")
map_out        <- file.path(dash_dir, "img", "maps")
for (d in c(data_dir, violin_rcp_out, violin_grp_out, map_out)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ==========================================================================
# 1. Species master list -- one row per species, read straight off the
#    per_species_sdm CSVs (the "base" file; species/code/group are identical
#    in the "anthro" file for the same species).
# ==========================================================================
cat("=== Building species list from per_species_sdm CSVs ===\n")

sdm_files <- list.files(sdm_dir, pattern = "_base_route_trends_sdm\\.csv$", full.names = TRUE)
species_meta <- map_dfr(sdm_files, function(f) {
  d <- read.csv(f, nrows = 1, stringsAsFactors = FALSE)
  tibble(code = d$species_code, name = d$species, group = d$group,
         slug = sub("_base_route_trends_sdm\\.csv$", "", basename(f)))
}) %>% arrange(name)

cat("Species found:", nrow(species_meta), "\n")
stopifnot(!any(duplicated(species_meta$code)))
stopifnot(!any(duplicated(species_meta$slug)))

# Map filenames use the species' common name with apostrophes stripped and
# spaces -> underscores (confirmed exact for all 469 species against the
# actual output/species_maps/ files before this script was written).
species_meta <- species_meta %>%
  mutate(map_slug = str_replace_all(str_remove_all(name, "'"), " ", "_"))

# ==========================================================================
# 2. Coverage check -- verify every expected source file exists before
#    copying anything, so a silent gap doesn't turn into a silently-broken
#    dashboard panel.
# ==========================================================================
cat("\n=== Checking source file coverage ===\n")

check_missing <- function(paths, label) {
  missing <- paths[!file.exists(paths)]
  cat(label, "-- missing:", length(missing), "of", length(paths), "\n")
  if (length(missing) > 0) print(head(missing, 10))
  missing
}

map_paths <- expand_grid(i = seq_len(nrow(species_meta)), model = models, scenario = scenarios) %>%
  mutate(path = file.path(maps_src_dir, model, scenario,
                          paste0(species_meta$map_slug[i], "_", model, "_", scenario, "_2010_2025.png")))
missing_maps <- check_missing(map_paths$path, "Range-shift maps")

rcp_violin_paths <- expand_grid(code = species_meta$code, model = models) %>%
  mutate(path = file.path(violin_dir, paste0("all_groups_", model, "_", code, "_trend_by_rcp_cld.png")))
missing_rcp_violins <- check_missing(rcp_violin_paths$path, "By-RCP violins")

grp_violin_paths <- expand_grid(code = species_meta$code, model = models) %>%
  mutate(path = file.path(violin_dir, paste0("all_groups_", model, "_", code, "_trend_by_group_cld.png")))
missing_grp_violins <- check_missing(grp_violin_paths$path, "By-group violins")
# Expected: Red-breasted Merganser (RBME) x2 models -- 0 routes fall in
# Contraction/Stable/Expansion for either scenario, so 4c never wrote this
# plot for it. Any OTHER species missing here is a real gap worth checking.
unexpected_missing_grp <- grp_violin_paths %>% filter(path %in% missing_grp_violins, code != "RBME")
if (nrow(unexpected_missing_grp) > 0) {
  cat("UNEXPECTED missing by-group violins (not RBME):\n")
  print(unexpected_missing_grp)
}

stopifnot("Missing map files -- run 3d_visualization_map.R first" = length(missing_maps) == 0)
stopifnot("Missing by-RCP violin files -- run 4c first" = length(missing_rcp_violins) == 0)
stopifnot("Unexpected missing by-group violin files" = nrow(unexpected_missing_grp) == 0)

# ==========================================================================
# 3. Table data from the two xlsx workbooks (RCP4.5 + RCP8.5 sheets each)
# ==========================================================================
cat("\n=== Reading xlsx tables ===\n")

read_grp_sheet <- function(model, scenario) {
  path <- file.path(files_dir, paste0("all_species_habitat_contraction_stable_expansion_2010_2025_", model, ".xlsx"))
  d <- openxlsx::read.xlsx(path, sheet = scenario_sheet[[scenario]])
  cat(" ", basename(path), scenario_sheet[[scenario]], "--", nrow(d), "rows\n")
  d
}

xlsx_raw <- expand_grid(model = models, scenario = scenarios) %>%
  mutate(data = map2(model, scenario, read_grp_sheet))

# tables[[code]][[model]][[scenario]] = list(Contraction=list(n,mean,sd)|NULL, Stable=..., Expansion=...)
# NOTE: openxlsx::read.xlsx() converts header cells like "Contraction n" into
# the column name "Contraction.n" (space -> dot) -- confirmed via
# colnames(read.xlsx(...)) before writing this, not guessed.
build_species_row <- function(row) {
  out <- list()
  for (g in grp_levels) {
    n <- row[[paste0(g, ".n")]]
    if (is.null(n) || is.na(n)) {
      out[[g]] <- NULL
    } else {
      out[[g]] <- list(n = unname(n),
                        mean = unname(row[[paste0(g, ".mean")]]),
                        sd = unname(row[[paste0(g, ".sd")]]))
    }
  }
  out
}

tables <- setNames(vector("list", nrow(species_meta)), species_meta$code)
for (code in species_meta$code) tables[[code]] <- setNames(vector("list", length(models)), models)

missing_from_xlsx <- character(0)
for (r in seq_len(nrow(xlsx_raw))) {
  model <- xlsx_raw$model[r]; scenario <- xlsx_raw$scenario[r]
  d <- xlsx_raw$data[[r]]
  d_by_name <- split(d, d$Species)
  for (i in seq_len(nrow(species_meta))) {
    code <- species_meta$code[i]; name <- species_meta$name[i]
    if (is.null(tables[[code]][[model]])) tables[[code]][[model]] <- list()
    if (!name %in% names(d_by_name)) {
      missing_from_xlsx <- c(missing_from_xlsx, paste(name, model, scenario))
      tables[[code]][[model]][[scenario]] <- list()  # empty -> "no data" in the UI
      next
    }
    tables[[code]][[model]][[scenario]] <- build_species_row(d_by_name[[name]][1, ])
  }
}
cat("\nSpecies/model/scenario combos with no xlsx row (expected: RBME x2 models x2 scenarios):",
    length(missing_from_xlsx), "\n")
if (length(missing_from_xlsx) > 0) print(missing_from_xlsx)

# ==========================================================================
# 4. Copy + downsize images to WebP
# ==========================================================================
cat("\n=== Converting images (this takes a few minutes) ===\n")

convert_one <- function(src, dst, width) {
  img <- magick::image_read(src)
  img <- magick::image_resize(img, paste0(width, "x"))
  magick::image_write(img, dst, format = "webp", quality = webp_quality)
}

t0 <- Sys.time()
n_done <- 0
report_every <- 200

for (i in seq_len(nrow(species_meta))) {
  code <- species_meta$code[i]
  for (model in models) {
    # By-RCP violin
    src <- file.path(violin_dir, paste0("all_groups_", model, "_", code, "_trend_by_rcp_cld.png"))
    convert_one(src, file.path(violin_rcp_out, paste0(code, "_", model, ".webp")), violin_width)
    n_done <- n_done + 1

    # By-group violin (skip for RBME, which has none)
    src <- file.path(violin_dir, paste0("all_groups_", model, "_", code, "_trend_by_group_cld.png"))
    if (file.exists(src)) {
      convert_one(src, file.path(violin_grp_out, paste0(code, "_", model, ".webp")), violin_width)
    }
    n_done <- n_done + 1

    for (scenario in scenarios) {
      src <- file.path(maps_src_dir, model, scenario,
                       paste0(species_meta$map_slug[i], "_", model, "_", scenario, "_2010_2025.png"))
      convert_one(src, file.path(map_out, paste0(code, "_", model, "_", scenario, ".webp")), map_width)
      n_done <- n_done + 1
    }
  }
  if (i %% report_every == 0) {
    cat(sprintf("  %d/%d species (%s elapsed)\n", i, nrow(species_meta),
               format(Sys.time() - t0, digits = 3)))
  }
}
cat(sprintf("Converted images for %d species in %s\n", nrow(species_meta), format(Sys.time() - t0, digits = 3)))

# ==========================================================================
# 5. Write manifest.js + tables.js
# ==========================================================================
cat("\n=== Writing manifest.js + tables.js ===\n")

manifest <- species_meta %>%
  mutate(has_group_violin = !(code == "RBME")) %>%   # only known gap
  select(code, name, group, has_group_violin) %>%
  arrange(name)

groups <- sort(unique(manifest$group))

manifest_json <- jsonlite::toJSON(
  list(species = manifest, groups = groups),
  auto_unbox = TRUE, pretty = FALSE, na = "null"
)
writeLines(paste0("const MANIFEST = ", manifest_json, ";"), file.path(data_dir, "manifest.js"))

tables_json <- jsonlite::toJSON(tables, auto_unbox = TRUE, pretty = FALSE, na = "null")
writeLines(paste0("const TABLES = ", tables_json, ";"), file.path(data_dir, "tables.js"))

cat("\n=== Dashboard data built ===\n")
cat("Species:", nrow(manifest), "\n")
cat("Open dashboard/index.html in a browser to view.\n")
