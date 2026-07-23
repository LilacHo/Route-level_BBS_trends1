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

