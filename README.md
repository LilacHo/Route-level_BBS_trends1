# Route-level_BBS_trends1

Route-level trend models for birds from the North American Breeding Bird Survey (BBS), using an intrinsic CAR (iCAR) spatial model fit in Stan via `bbsBayes2` and `cmdstanr`. Every species is fit two ways — a no-covariate `base` model and an `anthro` model with an anthropogenic land-use covariate on the slope — so route-level trends can be compared with and without that covariate, and each species' climate-suitability trajectory (SDM) can be related to its trend.

Adapted from Adam Smith's [Route-level_BBS_trends](https://github.com/AdamCSmithCWS/Route-level_BBS_trends).

## Model

For each species, route `r`, year `t`:

```
log(lambda[r,t]) = alpha[r] + beta[r]*t [+ gamma1 * anthro[r,t]]
```

`alpha[r]`/`beta[r]` (route-level intercept/slope) are given a non-centered iCAR (intrinsic CAR) spatial prior in Stan: the directly-sampled parameters are `alpha_raw[r]`/`beta_raw_space[r]`, and the reported `alpha[r]`/`beta[r]` are `generated quantities` — `alpha = sdalpha*alpha_raw + ALPHA`, `beta = sdbeta_space*beta_raw_space + BETA`. This distinction matters for convergence checking; see `route_converged` below.

- `models/slope_iCAR_route_NB_New.stan` — the `base` model (no covariate)
- `models/slope_iCAR_route_NB_New_covariate.stan` — adds `gamma1` on the `anthro` covariate

`base` and `anthro` are always fit on the IDENTICAL reduced dataset per species (only route-years where the anthro covariate is non-NA) — "same data for a fair comparison" is a rule followed throughout this pipeline, so any difference between the two models' trends reflects the covariate itself, not a difference in which rows were used.

## Pipeline

Run in order:

```
0_prepare_aou.R
  -> 1c_species_iCAR_covariates.R
       -> 1d_refit_nonconverged_species.R   (optional, only if some fits fail convergence)
  -> 2c_generate_route_trend_csvs_covariates.R
  -> 3c_add_SDM_covariates.R
  -> 4c_statistical_analysis_and_visualization_covariates.R
```

Every step downstream of 1c is a re-runnable post-processing step (safe to re-run at any point) except 1c/1d themselves, which are long model-fitting runs designed to be resumable across many interrupted sessions.

### 0_prepare_aou.R

Adds AOU numeric species codes to the target species list and checks each species against `bbsBayes2`'s species table.

- Reads: `data/spp_names_codes_group.csv`
- Writes: `data/spp_names_codes_group_aou.csv` (adds `spp.num`, `in_bbs`, `bbs_english`, etc. — `in_bbs == TRUE` is the species list every later script filters to)

### 1c_species_iCAR_covariates.R

Full production run: fits BOTH `base` and `anthro` for EVERY BBS species in the project's list (`in_bbs == TRUE`, currently ~607 species across all groups — aridlands, boreal/eastern/subtropical/western forests, grasslands, marshlands, coastal, waterbirds, generalists, arctic, urban/suburban), 2010-2025.

This replaced an earlier per-group, many-covariates design (a random species panel per land-cover group x several habitat/anthro/transition covariates). That pilot work found the `*_to_anthro` transition covariates were never identifiable, `anthro` was the one covariate that consistently produced a credible, well-powered effect, and group-specific habitat covariates weren't available consistently enough across all 12 groups the way `data/Anthro.csv` is. So this version drops the per-group machinery and just runs `anthro`, for every species, everywhere, flat (no per-group looping).

- Reads: `data/spp_names_codes_group_aou.csv`, `data/Anthro.csv` (one row per BBS route per year)
- Writes per species x tag: `output/rds/<species>_iCAR_<tag>_<firstYear>_<lastYear>_{stanfit,summ_fit}.rds`, `data/route_info/<species>_<tag>_..._route_info.rds`, `data/stan_data/<species>_<tag>_..._stan_data.RData`, plus a combined `output/diagnostics_covariates_all_species_anthro_<firstYear>_<lastYear>.csv`
- Resumable: `force_refit` defaults to `FALSE`, so a species/tag already fit is skipped — safe to stop and restart across many sessions for this ~600-species x 2-model run
- Always runs `helper/model_convergence.R`'s whole-model convergence table at the end (covers `base` too). Does NOT run `helper/gamma_lookup.R` automatically — run that separately when you want gamma1's own posterior summary/credibility.
- Shares its data-prep and fitting logic with `1d_refit_nonconverged_species.R` via `functions/covariate_model_fitting.R` (one implementation, not duplicated copies)

### 1d_refit_nonconverged_species.R

Targeted re-fit of ONLY the species/tag combinations that fail the whole-model convergence criterion (Rhat < 1.01 and bulk ESS > 400 across every parameter, `base` included) — not a blanket re-run of every species with longer iterations.

Steps, in order:

1. Reads every existing `*_summ_fit.rds` directly off disk (via `model_convergence_table()`) to find exactly which species/tag combos currently fail.
2. Re-fits only those, with increased iterations/adapt_delta (4000/4000, `adapt_delta = 0.9` vs. 1c's 2000/2000, 0.8) on the same reduced dataset 1c uses. Before overwriting, copies the current under-converged fit into parallel `output/rds_under-converged/`, `data/route_info_under-converged/`, `data/stan_data_under-converged/` backup folders (`archive_existing_fit()`), keeping only the most recent pre-refit snapshot.
3. Re-checks convergence after the standard refit.
4. If Western Meadowlark | anthro is still failing — a known marginal case where gamma1 itself sits just over threshold — gives it one more, more aggressive push (8000/8000, `adapt_delta = 0.95`) via `helper/western_meadowlark_refit.R`.
5. (Optional, commented out by default) route-level-exclusion drill-down: `helper/translate_offending_routes.R` then `helper/check_flagged_route_sparsity.R`, for whatever is still failing after step 4. Requires `helper/diagnose_nonconvergence.R`'s output to exist first.
6. Runs `helper/gamma_lookup.R`'s gamma1 lookup table LAST, once every model is set, so it reflects the post-refit fits.

Resumable like 1c — re-checks each species/tag's current status right before refitting, in case an earlier interrupted run of this script already fixed it. After running this, re-run `2c_generate_route_trend_csvs_covariates.R` so the improved fits propagate downstream.

### 2c_generate_route_trend_csvs_covariates.R

Post-processing only, no fitting. Combines per-route trend output across every species and both models (`base`, `anthro`) into a small number of CSVs, deriving `trend` (annual % change) and `rel_abundance` from each route's `alpha`/`beta` posterior means.

Also computes **`route_converged`**, a per-route boolean: TRUE only if THIS route's own `alpha_raw[r]` AND `beta_raw_space[r]` — the directly-sampled, non-centered parameters — both individually meet Rhat < 1.01 and bulk ESS > 400. This is deliberately based on the raw parameters, not the transformed `alpha[r]`/`beta[r]` generated quantities: when a species' variance hyperparameter (`sdalpha`/`sdbeta_space`) is small, a poorly-mixing raw parameter barely moves the transformed value, so the transformed parameter's own Rhat/ESS can look fine even while the underlying sampled parameter has real ICAR/funnel non-convergence. `route_converged` is independent of whole-model convergence (`model_converged`, from `helper/model_convergence.R`) — a route can be `TRUE` here even inside a model that fails whole-model convergence elsewhere, and it's this per-route flag, not the whole-model one, that drives downstream row exclusion.

- Reads: `output/rds/*_summ_fit.rds`, `data/route_info/*_route_info.rds` (both written by 1c/1d)
- Writes:
  - `output/species_routes_covariates/all_route_trends_all_species_anthro_<firstYear>_<lastYear>.csv` — one row per route/species/model: `species, species_code, group, model, route, latitude, longitude, alpha, beta, route_converged, trend, trend_lci, trend_uci, rel_abundance`
  - `output/species_routes_covariates/all_model_level_summary_all_species_anthro_<firstYear>_<lastYear>.csv` — one row per species/model: `species, species_code, group, model, n_routes, gamma1, gamma1_lci, gamma1_uci, gamma1_excludes_zero` (gamma1's 90% CI, `NA` for `base`)
  - `output/species_routes_covariates/per_species/<species>_<model_tag>_route_trends.csv` — same per-route columns, split one file per species x model tag, since `3c_add_SDM_covariates.R` requires exactly one `species_code` per input file

Safe to re-run at any point — species/tags 1c hasn't fit yet are skipped with a message, not an error.

### 3c_add_SDM_covariates.R

Adds `rcp45`/`rcp85` climate-scenario columns to each per-species-per-model route CSV from 2c, by extracting SDM classified-change raster values at each route's coordinates.

- Reads: `output/species_routes_covariates/per_species/*_route_trends.csv`
- Writes: `output/species_routes_covariates/per_species_sdm/*_route_trends_sdm.csv` (same columns, plus `rcp45`/`rcp85`)
- Raster path is resolved from each row's `group` column (written directly by 2c): `data/rcp{45,85}_<group>/<code>/<group>_<code>_breeding_2025_{45,85}_ENSEMBLE_classifiedchange.tif`
- A missing raster for one species/group is expected (not every group's rasters may exist yet) — that file is skipped with a warning, not treated as a fatal error
- Raster value legend (Bateman et al. 2020): 0 = never suitable, 1 = extirpation, 2 = worsening, 3 = slightly worsening, 4 = neutral, 5 = slightly improving, 6 = improving, 7 = colonization

### 4c_statistical_analysis_and_visualization_covariates.R

Statistical analysis and violin plots of route-level `trend` against SDM-derived range-change category, with compact letter displays (CLD) from pairwise Wilcoxon (BH-adjusted) tests.

- `model_tags` is hardcoded to `c("base", "anthro")` — Parts 1-4 below run once per tag in a loop (each writing its own tag-labeled files), so both models are always analyzed in one run.
- `bird_group` set to `NA` (default) pools every species/group across all Parts; set it to one real group to restrict.
- `require_route_converged` (default `TRUE`) drops rows whose `route_converged` isn't `TRUE` before any test or plot. Flipping this silently overwrites the previous run's output at the same filenames — rename/move it first if you want to keep both versions.

Parts, per model tag:

- **Part 1** — all 8 SDM categories (0-7), all species pooled: Kruskal-Wallis + pairwise Wilcoxon, violin plot with CLD.
- **Part 2** — same, per species.
- **Part 3** — grouped categories (Contraction = 1-3, Stable = 4, Expansion = 5-7), all species pooled.
- **Part 4** — same, per species.

After the loop:

- **Part 5** — base vs. anthro paired trend comparison. Since both models are fit on the identical reduced dataset per species, a route's `base` trend and `anthro` trend are a natural PAIRED comparison (same species, same route, same underlying counts). Uses a paired Wilcoxon signed-rank test overall and per species (BH-adjusted), writes a stats `.txt`, a per-species results `.csv`, and a histogram of the per-route trend difference.

- Reads: `output/species_routes_covariates/per_species_sdm/*_route_trends_sdm.csv`
- Writes: `output/species_routes_covariates/per_species_sdm_stats/` (stats `.txt`/`.csv`), `output/species_routes_covariates/per_species_sdm_plot/` (violin `.png`s, base-vs-anthro histogram)

## functions/ — shared library code

Pure function definitions, sourced (never run standalone) by the pipeline scripts above. If everything in here were deleted, the pipeline would break immediately — this is load-bearing, not optional.

- **`covariate_model_fitting.R`** — the single shared implementation of `species_to_f()`, `load_covariate()`, `prepare_species_data()` (BBS pull + covariate join + NA-drop + Voronoi spatial neighbours for one species), and `fit_one_covariate_model()` (fits `base` or a named covariate model for one species, saves stanfit/summary/route_info/stan_data, includes a chain-completion gate that surfaces a per-chain crash clearly instead of a cryptic downstream error). Used by both `1c` and `1d` so there is exactly one copy of this logic.
- **`posterior_summary_functions.R`** — general-purpose helpers (`posterior_samples()`, `posterior_sums()`) for tidying cmdstanr/rstan posterior draws into long/summary tables, carried over from the source project.
- **`neighbours_define_voronoi.R`** — builds the Voronoi-polygon spatial adjacency structure (clipped to each species' strata + concave hull of its routes) that feeds the iCAR prior's `node1`/`node2` edge list in Stan.

## helper/ — diagnostic tools

Standalone scripts for inspecting already-fitted output without re-fitting anything. Every setting is a function argument (nothing hardcoded to one bird group or species list), and each follows the same convention: sourcing the file with `<name>_skip_autorun <- TRUE` set beforehand loads just its function(s); sourcing it plain, or `Rscript helper/<name>.R`, also runs it standalone against this project's current defaults. Unlike `functions/`, nothing here is required for the main pipeline to run — this is the QA/inspection layer on top of it.

- **`model_convergence.R`** — whole-model convergence check (Rhat < 1.01, bulk ESS > 400 across every parameter) for every species x every model tag, `base` included. Reads `*_summ_fit.rds` directly off disk rather than 1c's own diagnostics CSV, since that CSV only logs species fit during one particular run session and this project's ~600-species run is designed to span many. This is `model_converged` — diagnostic only; it does not drive any downstream row exclusion (that's `route_converged`, computed in 2c).
- **`gamma_lookup.R`** — pulls gamma1's full posterior summary (mean/median/sd/q5/q95/rhat/ess) out of every species x model's `*_summ_fit.rds`, prints a credibility/convergence summary, writes `output/files/gamma_lookup_<label>_<firstYear>_<lastYear>.csv`.
- **`diagnose_nonconvergence.R`** — for every species/tag currently failing whole-model convergence, classifies which specific parameter(s) are driving the failure (`route_alpha`, `route_beta`, `gamma1`, or `hyperparameter_sd`/`other`) instead of just reporting the single worst Rhat/ESS. Distinguishes one sparsely-sampled route (often benign) from gamma1 or many routes being poorly mixed (more serious). Writes `output/files/nonconvergence_diagnosis_<firstYear>_<lastYear>.csv`.
- **`translate_offending_routes.R`** — translates `diagnose_nonconvergence.R`'s route-level offenders (internal `routeF` index) into real BBS route IDs and coordinates via each species/tag's own `route_info.rds`. The concrete follow-through on flagging/excluding specific routes rather than dropping an entire species. Writes `output/files/nonconvergence_flagged_routes_<firstYear>_<lastYear>.csv` (+ a per-species/tag summary).
- **`check_flagged_route_sparsity.R`** — for routes flagged by `translate_offending_routes.R`, tests whether they're more count-volatile (year-to-year CV) or more mismatched from their spatial neighbors' average count level than a species' other (converged) routes, via two-sample Wilcoxon rank-sum tests. Writes `output/files/flagged_route_sparsity_{check,detail}_<firstYear>_<lastYear>.csv`.
- **`western_meadowlark_refit.R`** — targeted extra-push refit (8000/8000 iterations, `adapt_delta = 0.95`) for Western Meadowlark | anthro specifically, a known marginal case where gamma1 itself remains just over threshold after the standard 1d refit. Called automatically from `1d_refit_nonconverged_species.R`'s Step 4; can also run standalone.
- **`route_info.R`** — backfills `data/route_info/*.rds` from an existing `stan_data.RData`, for any run where the route_info save didn't happen alongside the stan_data save.

## Convergence concepts (three distinct axes, deliberately not conflated)

- **`model_converged`** (`helper/model_convergence.R`) — whole-model, every parameter in a fit. Diagnostic only; does not itself drive any downstream row exclusion.
- **`route_converged`** (computed in `2c`) — per-route, based on that route's own `alpha_raw[r]`/`beta_raw_space[r]`. This IS what drives downstream filtering (`4c`'s `require_route_converged`).
- **gamma1's own credibility** — whether its 90% CI excludes zero (`gamma1_excludes_zero` in 2c's model-level CSV, and `helper/gamma_lookup.R`'s report). This is about statistical precision/significance, a separate question from MCMC convergence.

## Data layout

- `data/spp_names_codes_group_aou.csv` — species list with AOU codes and BBS-group membership (from `0_prepare_aou.R`)
- `data/Anthro.csv` — the anthro covariate, one row per BBS route per year
- `data/route_info/`, `data/stan_data/` — per-species/tag lightweight route lookup and full model input, written by 1c/1d
- `data/rcp45_<group>/`, `data/rcp85_<group>/` — SDM classified-change rasters, one folder per bird group
- `output/rds/` — per-species/tag `stanfit.rds`/`summ_fit.rds`
- `output/species_routes_covariates/` — 2c/3c/4c's combined and per-species CSVs, stats, and plots
- `models/` — the two Stan programs (`slope_iCAR_route_NB_New.stan`, `slope_iCAR_route_NB_New_covariate.stan`)
