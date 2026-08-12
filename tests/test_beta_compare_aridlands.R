## =============================================================================
## Test: compare route-level trend (beta[r]) between the "base" (no covariate)
## and "anthro" models for the 4 aridland species whose anthro gamma1 came
## back credible in the pilot panel (see output/files/aridlands_2010_2025.txt
## and the gamma_lookup_aridlands_2010_2025 write-up).
##
## Does NOT re-run 1c_species_iCAR_covariates_aridlands.R or refit anything —
## reads directly from the *_summ_fit.rds files 1c already wrote to
## output/rds/, via helper/beta_compare.R (which only needs dplyr + here, no
## cmdstanr/Stan). If those files aren't present locally (e.g. this is a
## fresh checkout, or the aridlands run hasn't been done on this machine
## yet), this will just report which files are missing rather than fitting
## anything.
##
## Purpose: beta with a covariate in the model is the trend LEFT OVER after
## the model reassigns whatever it can explain to the covariate (gamma1). If
## beta shrinks meaningfully from base -> anthro for a species, that's
## evidence anthro is explaining part of that species' raw decline; if beta
## barely moves, anthro isn't doing much of that species' story even though
## its gamma1 was credible in the panel-wide pilot.
## =============================================================================

library(here)

here::i_am("tests/test_beta_compare_aridlands.R")
source(here::here("helper", "beta_compare.R"))

firstYear <- 2010
lastYear  <- 2025
rds_dir   <- here::here("output", "rds")

## The 4 species with a credible aridland_anthro/anthro gamma1 in the pilot
## panel (Ladder-backed Woodpecker, Sage Thrasher, Bell's Vireo, Lark Sparrow
## -- see output/files/aridlands_2010_2025.txt, "anthro" section).
target_spp <- data.frame(
  Common.Name = c("Bell's Vireo",
                  "Sage Thrasher",
                  "Ladder-backed Woodpecker",
                  "Lark Sparrow")
)

beta_vs_anthro <- beta_compare_table(
  target_spp = target_spp,
  tag_a      = "base",
  tag_b      = "anthro",
  firstYear  = firstYear,
  lastYear   = lastYear,
  rds_dir    = rds_dir,
  label      = "aridlands_credible4",
  write_csv  = TRUE
)

## Uncomment to also check against aridland_anthro2 (the aridlands-specific
## anthro definition, distinct from the generic Anthro.csv used above) --
## same 4 species were credible under aridland_anthro2 in the pilot too:
# beta_vs_anthro2 <- beta_compare_table(
#   target_spp = target_spp,
#   tag_a      = "base",
#   tag_b      = "aridland_anthro2",
#   firstYear  = firstYear,
#   lastYear   = lastYear,
#   rds_dir    = rds_dir,
#   label      = "aridlands_credible4",
#   write_csv  = TRUE
# )

cat("\n=== Done ===\n")
if (!is.null(beta_vs_anthro)) {
  cat("Compared base vs anthro for", nrow(beta_vs_anthro), "of", nrow(target_spp), "species.\n")
} else {
  cat("No comparisons produced -- check that *_summ_fit.rds files for 'base' and\n",
      "'anthro' exist in", rds_dir, "for these 4 species.\n")
}
