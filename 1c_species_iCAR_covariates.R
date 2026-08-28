## =============================================================================
## Full production run: iCAR route-level slope model, base vs. anthro
## Time period: 2010-2025
##
## Fits TWO models for EVERY BBS species in this project's species list
## (data/spp_names_codes_group_aou.csv, in_bbs == TRUE -- all groups:
## aridlands, boreal_forests, eastern_forests, western_forests,
## subtropical_forests, grasslands, marshlands, coastal, waterbirds,
## generalists, arctic, urban_suburban -- currently 607 species):
##
##   base   -- no covariate    (models/slope_iCAR_route_NB_New.stan)
##   anthro -- data/Anthro.csv  (models/slope_iCAR_route_NB_New_covariate.stan)
##     log(lambda[r,t]) = alpha[r] + beta[r]*t + gamma1*Anthro[r,t]
##
## This replaces the earlier per-group, many-covariates design (random
## N-species panels x several habitat/anthro/transition covariates per
## forest group). That pilot work (see output/files/aridlands_2010_2025.txt
## and the gamma_lookup/beta_compare write-ups) found: (1) the "*_to_anthro"
## transition covariates never produced a credible gamma1 in any tested
## species -- too rare/sparse to be identifiable at this route/year
## resolution; (2) "anthro" is the one covariate that consistently produced
## a real, credible, well-powered effect; and (3) group-specific habitat
## covariates aren't available/consistent across every one of the 12 groups
## the way Anthro.csv is (it's already used project-wide). So this version
## drops the per-group covariate machinery entirely and just runs anthro,
## for everyone, everywhere.
##
## Both models are fit on the SAME reduced dataset per species (only
## route-years where Anthro is non-NA), same "same data for a fair
## comparison" principle used throughout this project.
##
## Covariate data: data/Anthro.csv, one row per BBS route per year
## (CountryNum/StateNum/Route/.../year/Anthro), loaded with load_covariate()
## the same way every other covariate in this project has been. Not
## rescaled (standing-stock proportion already spans most of [0,1]).
##
## Route-key format: new_data$route (from bbsBayes2::prepare_data()) is
## formatted "<StateNum>-<Route>" -- but ONLY when stratified by "bcr" (see
## every earlier version of this script for the full explanation). strat =
## "bcr" is used here for the same reason.
##
## Per user decision (carried over from every earlier version of this
## script): counts with a missing (NA) Anthro value are dropped from that
## species' stan_data, applied once up front, so both model tags see
## identical data.
##
## Species: EVERY species in the project's list, not a random panel --
## sorted alphabetically for a stable, resumable run order.
##
## force_refit defaults to FALSE here (unlike the earlier pilot versions of
## this script, which defaulted to TRUE) -- this is a full ~600-species
## production run, not a small re-run-every-time pilot, so by default this
## script SKIPS any species/tag that already has a saved summ_fit.rds and
## simply resumes where a previous, possibly-interrupted run left off. Set
## force_refit <- TRUE to force a clean re-fit of everything instead.
##
## This is a LONG run: ~600 species x 2 models x (BBS data pull + Voronoi
## neighbours + 4-chain/4000-iteration Stan fit) each. Expect this to take
## many hours to days depending on hardware -- that's expected, not a bug.
##
## This script fits the models and saves, per species and per model tag:
##   output/rds/<species>_iCAR_<tag>_<firstYear>_<lastYear>_stanfit.rds
##   output/rds/<species>_iCAR_<tag>_<firstYear>_<lastYear>_summ_fit.rds
##   data/route_info/<species>_<tag>_<firstYear>_<lastYear>_route_info.rds
##   data/stan_data/<species>_<tag>_<firstYear>_<lastYear>_stan_data.RData
## Plus one combined diagnostics CSV and (via helper/gamma_lookup.R) one
## gamma1 lookup table across all species.
## =============================================================================

library(bbsBayes2)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
library(spdep)
library(concaveman)
library(here)

here::i_am("1c_species_iCAR_covariates.R")

# prepare_species_data() and fit_one_covariate_model() -- and the
# neighbours_define_voronoi()/posterior_summary_functions.R sources they
# depend on -- now live in functions/covariate_model_fitting.R, shared with
# 1d_refit_nonconverged_species.R and tests/test_green_heron_debug.R so
# there is exactly one implementation instead of three hand-kept copies.
source(here::here("functions", "covariate_model_fitting.R"))

# Ensure cmdstanr writes output to CSV (default) and set output_dir for temp
# files. Each species/model_tag combo gets its own subfolder under here
# (created/deleted inside fit_one_covariate_model()) so raw per-chain CSVs
# don't pile up in %TEMP% across a run of ~600 species x 2 model tags.
cmdstanr_output_dir <- file.path(tempdir(), "cmdstan_output")
if (!dir.exists(cmdstanr_output_dir)) dir.create(cmdstanr_output_dir, recursive = TRUE)

# Settings ----------------------------------------------------------------
firstYear  <- 2010
lastYear   <- 2025
dt         <- lastYear - firstYear

strat <- "bcr"   # NOT "bbs_usgs" -- see earlier versions of this script for
                 # the full explanation.

# Re-fit control: FALSE by default for this full-species run -- see header
# comment. Set TRUE to force a clean re-fit of everything.
force_refit <- FALSE

model_tags <- c("base", "anthro")

# Output directories --------------------------------------------------------
output_dir <- here::here("output")
if (!dir.exists(output_dir)) dir.create(output_dir)
rds_dir <- here::here("output", "rds")
if (!dir.exists(rds_dir)) dir.create(rds_dir, recursive = TRUE)
if (!dir.exists(here::here("data"))) dir.create(here::here("data"))
if (!dir.exists(here::here("data", "maps"))) dir.create(here::here("data", "maps"), recursive = TRUE)
if (!dir.exists(here::here("data", "stan_data"))) dir.create(here::here("data", "stan_data"), recursive = TRUE)
route_info_dir <- here::here("data", "route_info")
if (!dir.exists(route_info_dir)) dir.create(route_info_dir, recursive = TRUE)

# load_covariate() now comes from functions/covariate_model_fitting.R --
# one row per route per year; builds a "<StateNum>-<Route>" key to match
# new_data$route below. out_name becomes both the joined column name in
# new_data AND the model_tag used to select it in fit_one_covariate_model().

# Load one covariate per the spec above, optionally rescaling by its own
# observed max. Kept as a thin wrapper (rather than calling load_covariate()
# directly) so this script stays easy to extend with another covariate
# later without restructuring the loop below.
load_group_covariate <- function(file, value_col, out_name, rescale) {
  lk <- load_covariate(file, value_col, out_name)
  if (rescale) {
    cscale <- max(lk[[out_name]], na.rm = TRUE)
    lk[[out_name]] <- lk[[out_name]] / cscale
    cat("  ", out_name, "covariate rescaled by its observed max (", round(cscale, 6), ") to span [0,1]\n")
  }
  cat("  ", out_name, "covariate:", nrow(lk), "route-year rows (",
      length(unique(lk$route_key)), "routes)\n")
  lk
}

# Covariate spec -- "anthro" doubles as both the model_tag AND the joined
# column name in new_data (see prepare_species_data()/
# fit_one_covariate_model() below, which key off this directly instead of a
# per-tag if/else chain). Anthro.csv is used project-wide, independent of
# bird group, so there is no per-group covariate config here anymore.
covariate_specs <- list(
  anthro = list(file = "Anthro.csv", value_col = "Anthro", rescale = FALSE)
)

# Species list -- EVERY species, not a random panel -----------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                   stringsAsFactors = FALSE)

target_spp <- spp_df %>%
  filter(in_bbs == TRUE) %>%
  distinct(Common.Name, Code, .keep_all = TRUE) %>%
  arrange(Common.Name)

cat("=== Full run: ", nrow(target_spp), " species x ", length(model_tags),
    " models (", paste(model_tags, collapse = ", "), ") ===\n", sep = "")
cat("Groups represented:", paste(sort(unique(target_spp$Group)), collapse = ", "), "\n")

# species_to_f() now comes from functions/covariate_model_fitting.R

# Compile the Stan models once (reused across every species) ---------------
model_base   <- cmdstan_model("models/slope_iCAR_route_NB_New.stan",
                              stanc_options = list("O1"))
model_single <- cmdstan_model("models/slope_iCAR_route_NB_New_covariate.stan",
                              stanc_options = list("O1"))

# Load the (single) covariate once, up front, for every species ------------
cat("Loading covariates...\n")
covariate_lookups <- list()
for (cn in setdiff(model_tags, "base")) {
  spec <- covariate_specs[[cn]]
  covariate_lookups[[cn]] <- load_group_covariate(spec$file, spec$value_col, cn, spec$rescale)
}

# prepare_species_data() and fit_one_covariate_model() now come from
# functions/covariate_model_fitting.R (see source() call near the top of
# this script) -- previously defined inline here, now shared with
# 1d_refit_nonconverged_species.R and tests/test_green_heron_debug.R.

# ==========================================================================
# Main loop: every species, both model tags. Flat -- no per-group looping
# anymore, since anthro applies the same way to every species regardless of
# its Group.
# ==========================================================================
# Hardcoded into this script's own flow (not optional/manual): both
# gamma_lookup_table() (gamma1's own posterior summary + credibility) and
# model_convergence_table() (whole-model Rhat/ESS pass-fail, base included)
# are ALWAYS run at the end of this script, right after the main loop below.
gamma_lookup_skip_autorun <- TRUE
source(here::here("helper", "gamma_lookup.R"))
model_convergence_skip_autorun <- TRUE
source(here::here("helper", "model_convergence.R"))

results_list <- list()
diagnostics_list <- list()

for (i in seq_len(nrow(target_spp))) {
  sp       <- target_spp$Common.Name[i]
  sp_f     <- species_to_f(sp)
  sp_code  <- target_spp$Code[i]
  sp_bbs   <- target_spp$bbs_english[i]
  sp_group <- target_spp$Group[i]

  cat("\n================================================================\n")
  cat("  [", i, "/", nrow(target_spp), "]", sp, " (", sp_group, ")\n")
  cat("================================================================\n")

  all_exist <- all(vapply(model_tags, function(tag) {
    out_base  <- paste0(sp_f, "_iCAR_", tag, "_", firstYear, "_", lastYear)
    summ_file <- file.path(rds_dir, paste0(out_base, "_summ_fit.rds"))
    file.exists(summ_file)
  }, logical(1)))

  if (all_exist && !force_refit) {
    cat("  Skipping (both", paste(model_tags, collapse = " and "), "already fitted for", sp, ")\n")
    next
  }

  prepped <- tryCatch(
    prepare_species_data(species = sp, species_bbs = sp_bbs, strat = strat,
                         firstYear = firstYear, lastYear = lastYear,
                         covariate_lookups = covariate_lookups),
    error = function(e) {
      message("  [ERROR] Data prep failed for ", sp, ": ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(prepped)) next

  for (tag in model_tags) {
    out_base       <- paste0(sp_f, "_iCAR_", tag, "_", firstYear, "_", lastYear)
    summ_file      <- file.path(rds_dir, paste0(out_base, "_summ_fit.rds"))
    stan_data_file <- here::here("data", "stan_data",
                                 paste0(sp_f, "_", tag, "_", firstYear, "_", lastYear, "_stan_data.RData"))
    if (file.exists(summ_file) && file.exists(stan_data_file) && !force_refit) {
      cat("  Skipping (already fitted):", basename(summ_file), "\n")
      next
    }

    diagnostics <- tryCatch(
      fit_one_covariate_model(species = sp, species_f = sp_f, model_tag = tag,
                              prepped = prepped, firstYear = firstYear, lastYear = lastYear,
                              model_base = model_base, model_single = model_single,
                              rds_dir = rds_dir, route_info_dir = route_info_dir,
                              cmdstanr_output_dir = cmdstanr_output_dir,
                              chains = 4, iter_warmup = 2000, iter_sampling = 2000,
                              adapt_delta = 0.8, max_treedepth = 10,
                              show_exceptions = FALSE),
      error = function(e) {
        message("  [ERROR] Model '", tag, "' failed for ", sp, ": ", conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(diagnostics)) next

    diagnostics$group <- sp_group
    results_list[[paste(sp, tag, sep = " | ")]] <- summ_file
    diagnostics_list[[paste(sp, tag, sep = " | ")]] <- diagnostics

    cat("  Done —", tag, "fit + diagnostics saved for", sp, "\n")
  }

  # Write/refresh the combined diagnostics CSV every 25 species, not just at
  # the very end -- with a ~600-species, multi-hour/day run, losing all
  # progress info to an interruption right before the final write would be
  # painful. Cheap to do (append-in-memory, rewrite whole file).
  if (length(diagnostics_list) > 0 && i %% 25 == 0) {
    diag_csv <- file.path(output_dir,
                          paste0("diagnostics_covariates_all_species_anthro_",
                                 firstYear, "_", lastYear, ".csv"))
    write.csv(bind_rows(diagnostics_list), diag_csv, row.names = FALSE)
    cat("  [checkpoint] diagnostics written to:", diag_csv, "(", i, "/", nrow(target_spp), "species processed )\n")
  }
}

# Final combined diagnostics CSV ---------------------------------------------
if (length(diagnostics_list) > 0) {
  diagnostics_all <- bind_rows(diagnostics_list)
  diag_csv <- file.path(output_dir,
                        paste0("diagnostics_covariates_all_species_anthro_",
                               firstYear, "_", lastYear, ".csv"))
  write.csv(diagnostics_all, diag_csv, row.names = FALSE)
  cat("\nDiagnostics written to:", diag_csv, "\n")
}

# gamma1 lookup table across every species (non-base tags only) --------------
gamma_lookup_table(target_spp = target_spp, model_tags = model_tags,
                   firstYear = firstYear, lastYear = lastYear,
                   rds_dir = rds_dir, bird_group = "all_species_anthro")

# Whole-model convergence table across every species AND every tag,
# including "base" -- Rhat < 1.01 and bulk ESS > 400 required across every
# parameter, not just gamma1 (see helper/model_convergence.R).
model_convergence_table(target_spp = target_spp, model_tags = model_tags,
                        firstYear = firstYear, lastYear = lastYear,
                        rds_dir = rds_dir, bird_group = "all_species_anthro")

cat("\n\n=== DONE ===\n")
cat("Species x model fits this run:", length(results_list), "\n")
cat("Fits saved in:", output_dir, "\n")
