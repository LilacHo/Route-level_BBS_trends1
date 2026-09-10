# Route-level_BBS_trends1

Route-level trend models for birds from the North American Breeding Bird Survey (BBS), using an intrinsic CAR (iCAR) spatial model fit in Stan via `bbsBayes2` and `cmdstanr`.

Adapted from Adam Smith's [Route-level_BBS_trends](https://github.com/AdamCSmithCWS/Route-level_BBS_trends). 

## Pipeline

Run in order:

1. `0_prepare_aou.R`
2. `1_species_iCAR_2010_2025.R`
3. `2_generate_route_trend_csvs.R`

### 0_prepare_aou.R

Adds AOU numeric species codes to the target species list and checks each species against `bbsBayes2`'s species table.

- Reads: `data/spp_names_codes_group.csv`
- Writes: `data/spp_names_codes_group_aou.csv`
- Original to this repo — no equivalent in the source project.

### 1_species_iCAR_2010_2025.R

Fits the iCAR route-level trend model, one species at a time, over 2010-2025. For each species it prepares BBS data with `bbsBayes2`, builds a route-neighbour (Voronoi) adjacency structure, fits the Stan model, saves the posterior summary and stan data, and runs a leave-future-out cross-validation on the final year.

- Reads: `data/spp_names_codes_group_aou.csv`
- Writes per species: `output/<species>_iCAR_New_2010_2025_stanfit.rds`, `..._summ_fit.rds`, `data/stan_data/<species>_2010_2025_stan_data.RData`, plus a combined `output/diagnostics_grasslands_2010_2025.csv`
- Adapted from `1alt_Species_data_prep_bbsBayes2.R` (bbsBayes2-based data prep) and `Fitting_new_iCAR_slope_model.R` (iCAR model) in the source repo, combined into a single fit script.
- Uses `functions/neighbours_define_voronoi.R` and `functions/posterior_summary_functions.R`, and the Stan models in `models/`, all carried over from the source repo.

### 2_generate_route_trend_csvs.R

Post-processing only — no model fitting. Reads the saved posterior summaries from step 1 and writes, per route, the raw `alpha` (log-scale intercept) and `beta` (log-scale slope) posterior means alongside their back-transformed, interpretable versions: `trend` (annual % change, with a 90% CI) and `rel_abundance`.

- Reads: `output/<species>_iCAR_New_2010_2025_summ_fit.rds`, `data/stan_data/<species>_2010_2025_stan_data.RData`
- Writes: `output/species_routes/<species>_route_trends.csv` with columns `species, species_code, route, latitude, longitude, alpha, beta, trend, trend_lci, trend_uci, rel_abundance`
- `trend = 100 * (exp(beta) - 1)`: converts the log-scale annual slope `beta` into an annual percent change (`trend_lci`/`trend_uci` are the same transform applied to the 5th/95th posterior percentiles).
- `rel_abundance = exp(alpha)`: converts the log-scale intercept `alpha` back to the count scale — the model's expected count for an average observer at that route in the fixed (mid-series) year.
- New script for this repo — the source project has no direct equivalent step.

## Covariate pipeline

Parallel to the base pipeline above, this repo also fits models with a land-cover covariate added, one bird group (grasslands, forests, aridlands, ...) at a time:

```
log(lambda[r,t]) = alpha[r] + beta[r]*t + gamma1 * covariate[r,t]
```

Run in order: `1c_species_iCAR_covariates*.R` -> `2c_generate_route_trend_csvs_covariates.R` -> `3c_add_SDM_covariates.R` -> `4c_statistical_analysis_and_visualization_covariates.R`.

### 1c_species_iCAR_covariates.R / 1c_species_iCAR_covariates_aridlands.R

Fits a "base" (no-covariate) model plus several single-covariate models per species, for one or more bird groups. Uses a generalized architecture: each `model_tag` doubles as the joined covariate column name in `new_data`, resolved through a `group_config` named list rather than hardcoded per group, so the same `prepare_species_data()`/`fit_one_covariate_model()` functions serve every group. Species panels are drawn as a reproducible random sample per group (`set.seed()`), not hand-picked.

- Reads: `data/spp_names_codes_group_aou.csv`, one covariate CSV per model_tag (route-year proportions, columns `CountryNum/StateNum/Route/.../year/<value col>`; `check.names = FALSE` preserves spaces in the header where present)
- Writes per species x model_tag: `output/rds/<species>_iCAR_<tag>_<firstYear>_<lastYear>_{stanfit,summ_fit}.rds`, `data/route_info/<species>_<tag>_<firstYear>_<lastYear>_route_info.rds`, `data/stan_data/<species>_<tag>_<firstYear>_<lastYear>_stan_data.RData`, plus a combined `output/diagnostics_covariates_<land_cover>_<firstYear>_<lastYear>.csv` and (via `helper/gamma_lookup.R`) `output/files/gamma_lookup_<land_cover>_<firstYear>_<lastYear>.csv`
- `1c_species_iCAR_covariates.R` currently covers four forest groups in one run (boreal_forests, eastern_forests, subtropical_forests, western_forests); `1c_species_iCAR_covariates_aridlands.R` is a single-group adaptation fitting 9 model tags for aridlands (base + 8 covariate variants: standing-stock habitat/anthro levels plus habitat-to-anthro transition covariates)
- Uses `models/slope_iCAR_route_NB_New.stan` (base) and `models/slope_iCAR_route_NB_New_covariate.stan` (single-covariate)

### 2c_generate_route_trend_csvs_covariates.R / 3c_add_SDM_covariates.R / 4c_statistical_analysis_and_visualization_covariates.R

Covariate-pipeline counterparts of `2_generate_route_trend_csvs.R` / `3_add_SDM.R` / `4_statistical_analysis_and_visualization.R` — same per-route trend derivation, climate-scenario extraction, and CLD-annotated statistical comparison, run against 1c's output (`output/rds/`, `data/route_info/`) instead of `1_species_iCAR_2010_2025.R`'s, with entirely separate output directories so the two pipelines never collide.

**Known gap:** these three scripts are still written for the original grasslands 4-tag design (`base`/`grassland_habitat`/`grassland_anthro`/`grassland_habitat_to_anthro`) and have not yet been generalized for the forest-group or aridlands model tags — running them against forest or aridlands output will need their `model_tags`/`land_cover` settings (and 2c's hardcoded 4-model assumption) updated first.

## Diagnostic helpers (helper/)

Small, soft-coded scripts for inspecting already-fitted output without re-running or refitting any model — every setting is a function argument, nothing is hardcoded to one bird group.

- **`gamma_lookup.R`** — pulls gamma1's full posterior summary (mean/median/sd/q5/q95/rhat/ess) out of every species x model_tag's `*_summ_fit.rds`, prints a credibility/convergence summary, writes `output/files/gamma_lookup_<land_cover>_<firstYear>_<lastYear>.csv`. Can be sourced standalone or called at the end of a 1c run.
- **`covariate_variance_decomp.R`** — for a route-year covariate CSV, splits its variance into between-route (absorbed by alpha[r], unusable for estimating trend) vs. within-route (the only part gamma1 can actually use) shares. Useful for spotting weakly-identified covariates before fitting — though note the within-route *share* (ICC) can be misleading when comparing covariates of very different overall scale; the absolute within-route variance is what predicts identifiability. Writes `output/files/covariate_variance_decomp_<label>_<firstYear>_<lastYear>.csv`.
- **`beta_compare.R`** — compares route-level trend (beta[r]) between two already-fitted model_tags for the same species (e.g. `"base"` vs. `"anthro"`), to see how much of the raw trend actually gets reassigned to a covariate once it's added, vs. how much is left over. Writes `output/files/beta_compare_<tag_a>_vs_<tag_b>_<label>_<firstYear>_<lastYear>.csv`.
- **`route_info.R`** — backfills `data/route_info/*.rds` from an existing `stan_data.RData`, for runs saved before 1c started writing `route_info.rds` directly.

## tests/

Exploratory/diagnostic scripts, not part of the main pipeline:

- `test_16covariates_pilot.R` — fits all 16 individual NLCD land-cover covariates simultaneously (one model per species, not one model per covariate) across several bird groups' pilot species panels.
- `test_covariate_scaling.R`, `test_grass_to_anthro_scaling.R`, `test_grass_to_anthro_scaling_multi_species.R` — compare covariate transformations (raw / standardized / sqrt-rescaled / within-between) for single- and multi-species fits.
- `test_beta_compare_aridlands.R` — runs `helper/beta_compare.R`'s base-vs-anthro comparison for the aridlands pilot species panel, reading already-saved fits with no refitting.

