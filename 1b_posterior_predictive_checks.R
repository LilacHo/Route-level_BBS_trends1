## =============================================================================
## Posterior predictive checks (PPC) for the iCAR route-level trend models
## Time period: 2010-2025
##
## Reads the per-species posterior draws + stan data saved by
## 1_species_iCAR_2010_2025.R:
##   output/rds/<species>_iCAR_NB_<firstYear>_<lastYear>_stanfit.rds
##   data/stan_data/<species>_<firstYear>_<lastYear>_stan_data.RData
##
## The fitted Stan model (models/slope_iCAR_route_NB_New.stan) has no
## y_rep/log_lik block, so replicated counts are simulated here in R from the
## saved posterior draws of the NB2 count model:
##
##   E_i     = beta[route_i]*(year_i - fixedyear) + alpha[route_i]
##              + sdobs*obs_raw[observer_i] + eta*firstyr_i
##   phi     = 1 / sqrt(sdnoise)
##   y_rep_i ~ NegBinomial(mu = exp(E_i), size = phi)
##
## which mirrors exactly the count ~ neg_binomial_2_log(E, phi) likelihood in
## the Stan model block.
##
## For each species, y_rep is compared to the observed counts (y) with:
##   - density overlay, ECDF overlay, and rootogram (bayesplot)
##   - Bayesian p-values for mean, sd, max, and proportion-of-zeros
##
## This is a read-only diagnostic step — it does not refit or alter any
## existing model output. It only adds new files:
##   output/ppc/<species>_ppc_dens.png
##   output/ppc/<species>_ppc_ecdf.png
##   output/ppc/<species>_ppc_rootogram.png
##   output/ppc/<species>_ppc_stats.png
##   output/ppc/ppc_bayes_pvalues_<land_cover>_<firstYear>_<lastYear>.csv
##
## Run this after 1_species_iCAR_2010_2025.R (needs its saved stanfit +
## stan_data output). Independent of 2/3, which only use the summary fit.
## =============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(here)

here::i_am("1b_posterior_predictive_checks.R")

# Settings (match 1_species_iCAR_2010_2025.R) ------------------------------
land_cover <- "grasslands"
firstYear  <- 2010
lastYear   <- 2025

# Number of posterior draws used to simulate y_rep (subsampled for speed)
n_draws <- 200
set.seed(123)

# Overwrite existing PPC plots/summary rows?
overwrite <- TRUE

# Directories
rds_dir <- here::here("output", "rds")   # all .rds fit output lives here
ppc_dir <- here::here("output", "ppc")
if (!dir.exists(ppc_dir)) dir.create(ppc_dir, recursive = TRUE)

# Target group species list -------------------------------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                   stringsAsFactors = FALSE)

target_spp <- spp_df %>%
  filter(Group == land_cover, in_bbs == TRUE) %>%
  distinct(Common.Name, Code, .keep_all = TRUE) %>%
  arrange(Common.Name)

cat("=== Posterior predictive checks ===\n")
cat("Group:", land_cover, " | Period:", firstYear, "-", lastYear, "\n")
cat("Species in group (n =", nrow(target_spp), ")\n")

# Helper: convert species name to file-safe format ---------------------------
species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

# Helper: sort the columns of a cmdstanr draws matrix by their array index,
# e.g. "beta[1]", "beta[2]", ... "beta[10]" -> ascending numeric order,
# so that column j lines up with route/observer index j.
sort_indexed_cols <- function(m) {
  idx <- as.integer(gsub(".*\\[|\\]", "", colnames(m)))
  m[, order(idx), drop = FALSE]
}

# Helper: simulate y_rep from a random subsample of posterior draws ---------
simulate_yrep <- function(stanfit, stan_data, n_draws = 200) {

  ndraws_total <- posterior::ndraws(stanfit$draws("BETA"))
  draw_ids <- sort(sample.int(ndraws_total, min(n_draws, ndraws_total)))

  beta_draws <- sort_indexed_cols(
    as.matrix(stanfit$draws("beta", format = "draws_matrix"))[draw_ids, , drop = FALSE])
  alpha_draws <- sort_indexed_cols(
    as.matrix(stanfit$draws("alpha", format = "draws_matrix"))[draw_ids, , drop = FALSE])
  obs_raw_draws <- sort_indexed_cols(
    as.matrix(stanfit$draws("obs_raw", format = "draws_matrix"))[draw_ids, , drop = FALSE])
  eta_draws     <- as.numeric(
    as.matrix(stanfit$draws("eta", format = "draws_matrix"))[draw_ids, ])
  sdobs_draws   <- as.numeric(
    as.matrix(stanfit$draws("sdobs", format = "draws_matrix"))[draw_ids, ])
  sdnoise_draws <- as.numeric(
    as.matrix(stanfit$draws("sdnoise", format = "draws_matrix"))[draw_ids, ])

  route     <- stan_data$route
  year      <- stan_data$year
  firstyr   <- stan_data$firstyr
  observer  <- stan_data$observer
  fixedyear <- stan_data$fixedyear
  ncounts   <- stan_data$ncounts

  yrep <- matrix(NA_real_, nrow = length(draw_ids), ncol = ncounts)

  for (s in seq_along(draw_ids)) {
    obs_s <- sdobs_draws[s] * obs_raw_draws[s, ]
    E_s <- beta_draws[s, route] * (year - fixedyear) +
      alpha_draws[s, route] +
      obs_s[observer] +
      eta_draws[s] * firstyr
    phi_s <- 1 / sqrt(sdnoise_draws[s])
    yrep[s, ] <- rnbinom(ncounts, size = phi_s, mu = exp(E_s))
  }

  list(yrep = yrep, draw_ids = draw_ids)
}

# Helper: Bayesian p-value for a test statistic T ---------------------------
# p = Pr(T(y_rep) >= T(y_obs)); values near 0 or 1 flag a misfitting T.
bayes_pval <- function(y, yrep, FUN) {
  t_obs <- FUN(y)
  t_rep <- apply(yrep, 1, FUN)
  mean(t_rep >= t_obs)
}

prop_zero <- function(x) mean(x == 0)

# ============================================================================
# Main loop: PPC plots + Bayesian p-values, one set per species
# ============================================================================
pval_list <- list()

for (i in seq_len(nrow(target_spp))) {
  sp      <- target_spp$Common.Name[i]
  sp_f    <- species_to_f(sp)
  sp_code <- target_spp$Code[i]

  cat("\n[", i, "/", nrow(target_spp), "]", sp, "\n")

  out_base       <- paste0(sp_f, "_iCAR_NB_", firstYear, "_", lastYear)
  stanfit_file   <- file.path(rds_dir, paste0(out_base, "_stanfit.rds"))
  stan_data_file <- here::here("data", "stan_data",
                               paste0(sp_f, "_", firstYear, "_", lastYear, "_stan_data.RData"))

  dens_png <- file.path(ppc_dir, paste0(sp_f, "_ppc_dens.png"))
  if (file.exists(dens_png) && !overwrite) {
    cat("  Skipping (already exists):", basename(dens_png), "\n")
    next
  }

  if (!file.exists(stanfit_file) || !file.exists(stan_data_file)) {
    cat("  No fitted output found — run 1_species_iCAR_2010_2025.R first. Skipping.\n")
    next
  }

  stanfit <- readRDS(stanfit_file)
  load(stan_data_file)   # loads stan_data, new_data, route_map, ...

  y <- stan_data$count

  sim <- tryCatch(simulate_yrep(stanfit, stan_data, n_draws = n_draws),
                   error = function(e) {
                     message("  [ERROR] PPC simulation failed for ", sp, ": ", conditionMessage(e))
                     return(NULL)
                   })
  if (is.null(sim)) next
  yrep <- sim$yrep

  # Plots (bayesplot) ---------------------------------------------------------
  yrep_plot <- yrep[seq_len(min(50, nrow(yrep))), , drop = FALSE]

  p_dens <- bayesplot::ppc_dens_overlay(y = log1p(y), yrep = log1p(yrep_plot)) +
    ggplot2::labs(title = paste(sp, "— PPC density overlay (log1p count)"))
  ggplot2::ggsave(dens_png, p_dens, width = 7, height = 5, dpi = 150)

  p_ecdf <- bayesplot::ppc_ecdf_overlay(y = y, yrep = yrep_plot) +
    ggplot2::labs(title = paste(sp, "— PPC ECDF overlay"))
  ggplot2::ggsave(file.path(ppc_dir, paste0(sp_f, "_ppc_ecdf.png")),
                   p_ecdf, width = 7, height = 5, dpi = 150)

  p_root <- tryCatch(
    bayesplot::ppc_rootogram(y = y, yrep = yrep, style = "hanging") +
      ggplot2::labs(title = paste(sp, "— PPC rootogram")),
    error = function(e) NULL)
  if (!is.null(p_root)) {
    ggplot2::ggsave(file.path(ppc_dir, paste0(sp_f, "_ppc_rootogram.png")),
                     p_root, width = 7, height = 5, dpi = 150)
  }

  p_stats <- bayesplot::bayesplot_grid(
    bayesplot::ppc_stat(y, yrep, stat = "mean") + ggplot2::labs(subtitle = "mean"),
    bayesplot::ppc_stat(y, yrep, stat = "sd")   + ggplot2::labs(subtitle = "sd"),
    bayesplot::ppc_stat(y, yrep, stat = "max")  + ggplot2::labs(subtitle = "max"),
    bayesplot::ppc_stat(y, yrep, stat = prop_zero) + ggplot2::labs(subtitle = "prop. zero"),
    titles = rep("", 4)
  )
  ggplot2::ggsave(file.path(ppc_dir, paste0(sp_f, "_ppc_stats.png")),
                   p_stats, width = 9, height = 7, dpi = 150)

  # Bayesian p-values -----------------------------------------------------
  pvals <- data.frame(
    species      = sp,
    species_code = sp_code,
    p_mean       = round(bayes_pval(y, yrep, mean), 3),
    p_sd         = round(bayes_pval(y, yrep, sd), 3),
    p_max        = round(bayes_pval(y, yrep, max), 3),
    p_prop_zero  = round(bayes_pval(y, yrep, prop_zero), 3)
  )
  pval_list[[sp]] <- pvals

  cat("  Bayesian p-values — mean:", pvals$p_mean, " sd:", pvals$p_sd,
      " max:", pvals$p_max, " prop.zero:", pvals$p_prop_zero, "\n")
  cat("  PPC plots ->", ppc_dir, "\n")
}

# Combined Bayesian p-value summary CSV --------------------------------------
if (length(pval_list) > 0) {
  pvals_all <- bind_rows(pval_list)
  pval_csv <- file.path(ppc_dir,
                        paste0("ppc_bayes_pvalues_", land_cover, "_",
                               firstYear, "_", lastYear, ".csv"))
  write.csv(pvals_all, pval_csv, row.names = FALSE)
  cat("\nBayesian p-values written to:", pval_csv, "\n")

  # Flag species where any test statistic is extreme (p < 0.05 or p > 0.95),
  # i.e. the observed data falls in the tail of its own posterior predictive
  # distribution for that statistic — a sign of possible model misfit.
  flagged <- pvals_all %>%
    filter(if_any(starts_with("p_"), ~ .x < 0.05 | .x > 0.95))
  if (nrow(flagged) > 0) {
    cat("\nSpecies with extreme Bayesian p-values (possible misfit):\n")
    print(flagged)
  } else {
    cat("\nNo species flagged with extreme Bayesian p-values.\n")
  }
}

cat("\n=== Summary ===\n")
cat("Species checked this run:", length(pval_list), "\n")
cat("PPC outputs in:", ppc_dir, "\n")
