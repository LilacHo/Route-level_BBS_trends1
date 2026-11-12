# 4h_sdm_tracking_stage1_draws.R
#
# Stage 1 of the Bayesian "SDM tracking" analysis. Run this where the Stan
# posterior draws (*_stanfit.rds from 1c/1d) are stored; stage 2
# (4i_sdm_tracking_stage2_bayes.R) only needs this script's small CSV output.
#
# Question: does SDM-projected habitat-suitability change track observed
# population trends, and for which species / groups does it not?
#
# For each species, in EVERY retained posterior draw of the route-level slopes
# beta[r] (anthro iCAR model), this script:
#   1. converts each route's slope to a trend, 100 * (exp(beta) - 1) %/yr;
#   2. regresses route trend on the route's SDM class score (ordinary least
#      squares across the species' routes, each route weighted equally) --
#      the TRACKING SLOPE, in %/yr per class step. Positive = trends are better
#      where the SDM projects improving suitability (the species tracks it);
#   3. also computes the slope separately on the worsening side (classes 1-4)
#      and the improving side (classes 4-7), and the mean trend in Contraction
#      (1-3), Stable (4) and Expansion (5-7) routes, and Contraction - Expansion.
# Repeating this over the draws gives each species a full posterior for each
# quantity. Because the draws are joint across routes, these posteriors carry
# both the route-level trend uncertainty and the spatial (iCAR) correlation
# between neighbouring routes -- the two things the glmmTMB analysis (4f) on
# posterior means could not.
#
# SDM class score (class_score below): ordinal steps centred on Neutral,
#   1 extirpation -3 | 2 worsening -2 | 3 slightly worsening -1 | 4 neutral 0 |
#   5 slightly improving +1 | 6 improving +2 | 7 colonization +3.
# The percent-change bands behind the classes are open-ended at the top
# (improving = 50%+, colonization = newly suitable), so equal steps are used
# rather than invented midpoints. Change class_score to try another scoring.
# Class 0 (never suitable) and routes without a class are excluded, as in 4f.
#
# Species need routes in at least 2 distinct classes for the overall slope,
# and at least 2 distinct classes within 1-4 (or 4-7) for each side slope.
#
# Reads:  <rds_dir>/<species>_iCAR_<tag>_<firstYear>_<lastYear>_stanfit.rds   (cmdstanr save_object)
#         <route_info_dir>/<species>_<tag>_<firstYear>_<lastYear>_route_info.rds (route <-> routeF)
#         output/species_routes_covariates/per_species_sdm/<species>_<tag>_route_trends_sdm.csv
#           (route, SDM class, route_converged; from 3c)
# Writes: output/files/sdm_tracking/sdm_tracking_stage1_species_<tag>_<scenario>_<firstYear>_<lastYear>.csv
#           one row per species: status, route counts per class, and for each
#           quantity its posterior mean, sd, 2.5% / 97.5% quantiles and P(> 0)
#         output/files/sdm_tracking/sdm_tracking_stage1_draws_<tag>_<scenario>_<firstYear>_<lastYear>.rds
#           (optional, save_draws) thinned draws per species, for later checks
#         output/files/sdm_tracking/sdm_tracking_stage1_log_<tag>_<scenario>_<firstYear>_<lastYear>.txt
#
# Resumable: species already finished (status "ok") in an existing output CSV
# are skipped, so the script can be stopped and restarted.
#
# Needs: tidyverse, posterior, cmdstanr (to read the saved CmdStanMCMC objects).

library(here)
library(tidyverse)
library(posterior)
library(cmdstanr)

here::i_am("4h_sdm_tracking_stage1_draws.R")

# Settings -------------------------------------------------------------------
model_tag  <- "anthro"
scenario   <- "rcp45"
firstYear  <- 2010
lastYear   <- 2025
rds_dir        <- here::here("output", "rds")          # matches 1c_species_iCAR_covariates.R (*_stanfit.rds)
route_info_dir <- here::here("data", "route_info")     # matches 1c_species_iCAR_covariates.R (*_route_info.rds)
sdm_dir        <- here::here("output", "species_routes_covariates", "per_species_sdm")
out_dir        <- here::here("output", "files", "sdm_tracking")

n_draws    <- 2000      # posterior draws kept per species (random subset of the 8000)
save_draws <- TRUE      # also save a thinned set of per-species draws (n_draws_saved)
n_draws_saved <- 500
seed       <- 20261014
class_score <- c(`1` = -3, `2` = -2, `3` = -1, `4` = 0, `5` = 1, `6` = 2, `7` = 3)
only_species <- NULL    # e.g. c("Song Sparrow", "Wild Turkey") for a quick test run
# SDM classes to drop before computing slopes (sensitivity analysis). integer(0) = main
# analysis (all classes 1-7); c(7) = drop colonization. A non-empty setting adds a tag such
# as "_no7" to every output file name, so it never overwrites the main results.
# Must match the same setting in 4i_sdm_tracking_stage2_bayes.R.
exclude_classes <- c(7)

run_tag  <- if (length(exclude_classes) > 0) paste0("_no", paste(sort(exclude_classes), collapse = "")) else ""
suffix   <- paste0(paste(model_tag, scenario, firstYear, lastYear, sep = "_"), run_tag)
out_csv  <- file.path(out_dir, paste0("sdm_tracking_stage1_species_", suffix, ".csv"))
out_rds  <- file.path(out_dir, paste0("sdm_tracking_stage1_draws_", suffix, ".rds"))
log_file <- file.path(out_dir, paste0("sdm_tracking_stage1_log_", suffix, ".txt"))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

species_to_f <- function(sp) gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", paste0(...))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

# OLS slope of each draw's route trends on x: one value per draw
slope_draws <- function(Y, x) {
  xc <- x - mean(x)
  as.vector(Y %*% xc) / sum(xc^2)
}
summ <- function(v, nm) {
  if (is.null(v) || all(is.na(v))) {
    return(setNames(as.list(rep(NA_real_, 5)), paste0(nm, c("_mean", "_sd", "_q025", "_q975", "_p_pos"))))
  }
  setNames(list(mean(v), sd(v), unname(quantile(v, 0.025)), unname(quantile(v, 0.975)), mean(v > 0)),
           paste0(nm, c("_mean", "_sd", "_q025", "_q975", "_p_pos")))
}

# One species ----------------------------------------------------------------
process_species <- function(sdm_file) {
  sdm_all <- read.csv(sdm_file, stringsAsFactors = FALSE)
  sp   <- sdm_all$species[1]
  code <- sdm_all$species_code[1]
  base <- tibble(species = sp, species_code = code)

  sdm <- sdm_all %>%
    filter(route_converged %in% TRUE, .data[[scenario]] %in% setdiff(1:7, exclude_classes)) %>%
    transmute(route, cls = as.integer(.data[[scenario]]))
  counts <- setNames(as.list(tabulate(sdm$cls, nbins = 7)), paste0("n_class", 1:7))
  base <- bind_cols(base, as_tibble(counts), tibble(n_routes = nrow(sdm), n_classes = n_distinct(sdm$cls)))
  if (n_distinct(sdm$cls) < 2) return(list(row = mutate(base, status = "fewer than 2 SDM classes"), draws = NULL))

  sp_f <- species_to_f(sp)
  fit_file <- file.path(rds_dir, paste0(sp_f, "_iCAR_", model_tag, "_", firstYear, "_", lastYear, "_stanfit.rds"))
  ri_file  <- file.path(route_info_dir, paste0(sp_f, "_", model_tag, "_", firstYear, "_", lastYear, "_route_info.rds"))
  if (!file.exists(fit_file)) return(list(row = mutate(base, status = "stanfit file not found"), draws = NULL))
  if (!file.exists(ri_file))  return(list(row = mutate(base, status = "route_info file not found"), draws = NULL))

  ri <- readRDS(ri_file)
  sdm <- sdm %>% inner_join(ri %>% distinct(route, routeF), by = "route")
  if (n_distinct(sdm$cls) < 2) return(list(row = mutate(base, status = "fewer than 2 SDM classes after route match"), draws = NULL))

  fit <- readRDS(fit_file)
  dm  <- fit$draws(variables = "beta", format = "draws_matrix")
  rm(fit); invisible(gc())
  set.seed(seed)
  keep <- if (nrow(dm) > n_draws) sort(sample(nrow(dm), n_draws)) else seq_len(nrow(dm))
  cols <- paste0("beta[", sdm$routeF, "]")
  miss <- setdiff(cols, colnames(dm))
  if (length(miss) > 0) return(list(row = mutate(base, status = paste("beta columns missing:", length(miss))), draws = NULL))
  B <- unclass(dm)[keep, cols, drop = FALSE]
  rm(dm); invisible(gc())

  Y <- 100 * (exp(B) - 1)                         # draws x routes, %/yr
  x <- unname(class_score[as.character(sdm$cls)])

  s_all <- slope_draws(Y, x)
  low <- sdm$cls <= 4; up <- sdm$cls >= 4
  s_low <- if (n_distinct(x[low]) >= 2) slope_draws(Y[, low, drop = FALSE], x[low]) else NULL
  s_up  <- if (n_distinct(x[up])  >= 2) slope_draws(Y[, up,  drop = FALSE], x[up])  else NULL
  mean_in <- function(idx) if (any(idx)) rowMeans(Y[, idx, drop = FALSE]) else NULL
  m_con <- mean_in(sdm$cls <= 3); m_sta <- mean_in(sdm$cls == 4); m_exp <- mean_in(sdm$cls >= 5)
  ce <- if (!is.null(m_con) && !is.null(m_exp)) m_con - m_exp else NULL

  row <- bind_cols(base, as_tibble(c(
    summ(s_all, "slope"), summ(s_low, "slope_low"), summ(s_up, "slope_up"),
    summ(m_con, "mean_con"), summ(m_sta, "mean_stable"), summ(m_exp, "mean_exp"), summ(ce, "ce"),
    list(cor_low_up = if (!is.null(s_low) && !is.null(s_up)) cor(s_low, s_up) else NA_real_,
         n_draws_used = length(keep))
  ))) %>% mutate(status = "ok")

  dr <- NULL
  if (save_draws) {
    k <- seq_len(min(n_draws_saved, length(s_all)))
    dr <- list(slope = s_all[k], slope_low = s_low[k], slope_up = s_up[k], ce = ce[k])
  }
  list(row = row, draws = dr)
}

# Run ------------------------------------------------------------------------
sdm_files <- list.files(sdm_dir, pattern = paste0("_", model_tag, "_route_trends_sdm\\.csv$"), full.names = TRUE)
if (length(sdm_files) == 0) stop("No *_", model_tag, "_route_trends_sdm.csv files in ", sdm_dir)
if (!is.null(only_species)) {
  sdm_files <- sdm_files[basename(sdm_files) %in% paste0(species_to_f(only_species), "_", model_tag, "_route_trends_sdm.csv")]
}

done <- tibble(species = character(), status = character())
results <- list(); draws_out <- list()
if (file.exists(out_csv)) {
  prev <- read.csv(out_csv, stringsAsFactors = FALSE)
  results <- split(prev, seq_len(nrow(prev))) %>% map(as_tibble)
  done <- prev %>% select(species, status)
  if (save_draws && file.exists(out_rds)) draws_out <- readRDS(out_rds)
  log_msg("Resuming: ", sum(done$status == "ok"), " species already done")
}

log_msg("Stage 1 start: ", length(sdm_files), " species files; draws from ", rds_dir)
for (i in seq_along(sdm_files)) {
  f <- sdm_files[i]
  sp <- read.csv(f, nrows = 1, stringsAsFactors = FALSE)$species
  if (sp %in% done$species[done$status == "ok"]) next
  res <- tryCatch(process_species(f), error = function(e) {
    list(row = tibble(species = sp, status = paste("error:", conditionMessage(e))), draws = NULL)
  })
  results <- results[!vapply(results, function(r) identical(r$species[1], sp), logical(1))]
  results[[length(results) + 1]] <- res$row
  if (!is.null(res$draws)) draws_out[[res$row$species_code[1]]] <- res$draws
  log_msg(sprintf("[%d/%d] %s: %s%s", i, length(sdm_files), sp, res$row$status[1],
                  if (res$row$status[1] == "ok") sprintf(" (slope %.2f, sd %.2f)", res$row$slope_mean, res$row$slope_sd) else ""))
  if (i %% 25 == 0 || i == length(sdm_files)) {          # checkpoint
    write.csv(bind_rows(results), out_csv, row.names = FALSE)
    if (save_draws) saveRDS(draws_out, out_rds)
  }
}
out <- bind_rows(results)
write.csv(out, out_csv, row.names = FALSE)
if (save_draws) saveRDS(draws_out, out_rds)

# Short summary for the log ------------------------------------------------------
ok <- out %>% filter(status == "ok")
log_msg("Stage 1 done. Species ok: ", nrow(ok), " of ", nrow(out))
print(table(out$status))
log_msg("Status counts: ", paste(names(table(out$status)), table(out$status), sep = " = ", collapse = "; "))
if (nrow(ok) > 0) {
  log_msg(sprintf("Tracking slope across species (raw, unpooled): median %.3f; %d of %d species with P(slope > 0) > 0.9; %d with < 0.1",
                  median(ok$slope_mean), sum(ok$slope_p_pos > 0.9), nrow(ok), sum(ok$slope_p_pos < 0.1)))
}
cat("\nWrote:\n ", out_csv, "\n ", if (save_draws) out_rds else "", "\n ", log_file, "\n")
