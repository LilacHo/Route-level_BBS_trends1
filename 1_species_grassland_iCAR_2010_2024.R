## =============================================================================
## Grassland birds: New iCAR route-level trend model
## Time period: 2010-2024
## Loops through all grassland species in spp_names_codes_group.csv
## Checks each species exists in bbsBayes2 before fitting
## Based on Fitting_new_iCAR_slope_model.R (sum_to_zero_vector parameterization)
## =============================================================================

library(bbsBayes2)
library(tidyverse)
library(cmdstanr)
library(sf)
library(spdep)
library(concaveman)
library(patchwork)

source("functions/neighbours_define_voronoi.R")
source("functions/posterior_summary_functions.R")

output_dir <- "output"
if (!dir.exists(output_dir)) dir.create(output_dir)
if (!dir.exists("data")) dir.create("data")
if (!dir.exists("data/maps")) dir.create("data/maps", recursive = TRUE)
if (!dir.exists("figures")) dir.create("figures")

strat <- "bbs_usgs"
firstYear <- 2010
lastYear <- 2024

## =============================================================================
## Load species list and filter to grassland birds
## =============================================================================

spp_df <- read.csv("data/spp_names_codes_group.csv", stringsAsFactors = FALSE)
grassland_spp <- spp_df %>%
  filter(Group == "grasslands") %>%
  distinct(Common.Name, .keep_all = TRUE) %>%
  arrange(Common.Name)

cat("=== Grassland species iCAR model fitting ===\n")
cat("Period:", firstYear, "-", lastYear, "\n")
cat("Grassland species found in CSV:", nrow(grassland_spp), "\n\n")

## =============================================================================
## Load BBS data once and get available species list
## =============================================================================

strat_data <- load_bbs_data()
bbs_species <- strat_data$species %>%
  filter(unid_combined == TRUE)

# Check which grassland species exist in bbsBayes2
grassland_spp <- grassland_spp %>%
  mutate(in_bbs = Common.Name %in% bbs_species$english)

cat("Species available in bbsBayes2:",
    sum(grassland_spp$in_bbs), "/", nrow(grassland_spp), "\n")

if (any(!grassland_spp$in_bbs)) {
  cat("\nSpecies NOT found in bbsBayes2 (will be skipped):\n")
  cat(paste(" -", grassland_spp$Common.Name[!grassland_spp$in_bbs]), sep = "\n")
  cat("\n")
}

species_list <- grassland_spp$Common.Name[grassland_spp$in_bbs]

# Compile the Stan model once (reused for all species)
mod.file <- "models/slope_iCAR_route_NB_New.stan"
slope_model <- cmdstan_model(mod.file, stanc_options = list("O1"))

cv_mod.file <- "models/slope_iCAR_route_NB_cv.stan"
cv_model <- cmdstan_model(cv_mod.file, stanc_options = list("O1"))

# Track results across species
all_results <- NULL

## =============================================================================
## SPECIES LOOP
## =============================================================================

for (species in species_list) {

  species_f <- gsub(gsub(species, pattern = " ", replacement = "_", fixed = TRUE),
                    pattern = "'", replacement = "", fixed = TRUE)

  cat("\n================================================================\n")
  cat("  ", species, " (", which(species_list == species), "/",
      length(species_list), ")\n")
  cat("================================================================\n")

  # Check if already fitted — skip if output exists
  out_summ_file <- paste0(output_dir, "/", species_f, "_iCAR_New_",
                          firstYear, "_", lastYear, "_summ_fit.rds")
  if (file.exists(out_summ_file)) {
    cat("  Already fitted — skipping.\n")
    next
  }

  ## --- Data preparation ---
  tryCatch({

    data_pkg <- bbsBayes2::stratify(by = strat, species = species,
                                     use_map = FALSE) %>%
      bbsBayes2::prepare_data(min_year = firstYear,
                              max_year = lastYear,
                              min_n_routes = 1,
                              min_max_route_years = 1)

    raw_data <- data_pkg[["raw_data"]]
    strata_map <- load_map(strat)

    route_map1 <- raw_data %>%
      select(route, strata_name, latitude, longitude) %>%
      distinct()
    route_map1 <- st_as_sf(route_map1, coords = c("longitude", "latitude"))
    st_crs(route_map1) <- 4326
    route_map1 <- st_transform(route_map1, crs = st_crs(strata_map))

    strata_map_buf <- strata_map %>%
      filter(strata_name %in% route_map1$strata_name) %>%
      summarise() %>%
      st_buffer(., 10000)

    realized_routes <- route_map1 %>%
      st_join(., strata_map_buf, join = st_within, left = FALSE)

    new_data <- data.frame(
      strat_name = raw_data$strata_name,
      strat      = raw_data$strata,
      route      = raw_data$route,
      latitude   = raw_data$latitude,
      longitude  = raw_data$longitude,
      count      = raw_data$count,
      year       = raw_data$year_num,
      firstyr    = raw_data$first_year,
      ObsN       = raw_data$observer,
      r_year     = raw_data$year
    ) %>%
      filter(route %in% realized_routes$route)

    if (nrow(new_data) < 50) {
      cat("  Too few observations (", nrow(new_data), ") — skipping.\n")
      next
    }

    realized_strata_map <- strata_map %>%
      filter(strata_name %in% unique(new_data$strat_name))

    new_data$routeF <- as.integer(factor(new_data$route))

    route_map <- unique(data.frame(
      route     = new_data$route,
      routeF    = new_data$routeF,
      strat     = new_data$strat_name,
      latitude  = new_data$latitude,
      longitude = new_data$longitude
    ))

    dups <- which(duplicated(route_map[, c("latitude", "longitude")]))
    while (length(dups) > 0) {
      route_map[dups, "latitude"]  <- route_map[dups, "latitude"]  + 0.01
      route_map[dups, "longitude"] <- route_map[dups, "longitude"] + 0.01
      dups <- which(duplicated(route_map[, c("latitude", "longitude")]))
    }

    route_map <- st_as_sf(route_map, coords = c("longitude", "latitude"))
    st_crs(route_map) <- 4326
    route_map <- st_transform(route_map, crs = st_crs(strata_map))

    car_stan_dat <- neighbours_define_voronoi(
      real_point_map = route_map,
      species = species,
      strat_indicator = "routeF",
      strata_map = realized_strata_map,
      concavity = 1,
      save_plot_data = TRUE,
      plot_dir = "data"
    )

    cat("  Routes:", max(new_data$routeF),
        " | Obs:", nrow(new_data),
        " | Edges:", car_stan_dat$N_edges, "\n")

    ## --- Fit iCAR model ---
    tmp <- data.frame(route = new_data$routeF, count = new_data$count) %>%
      group_by(route) %>%
      summarise(mean_count = mean(log(count + 1)))
    sd_alpha_prior <- sd(tmp$mean_count)

    stan_data <- list(
      count      = new_data$count,
      ncounts    = nrow(new_data),
      strat      = new_data$strat,
      route      = new_data$routeF,
      year       = new_data$year,
      firstyr    = new_data$firstyr,
      fixedyear  = floor(mean(1:max(new_data$year))),
      nyears     = max(new_data$year),
      observer   = as.integer(factor(new_data$ObsN)),
      nobservers = length(unique(new_data$ObsN)),
      N_edges    = car_stan_dat$N_edges,
      node1      = car_stan_dat$node1,
      node2      = car_stan_dat$node2,
      nroutes    = max(new_data$routeF),
      sd_alpha_prior = sd_alpha_prior
    )

    if (car_stan_dat$N != stan_data[["nroutes"]]) {
      cat("  Route mismatch in adjacency matrix — skipping.\n")
      next
    }

    stanfit <- slope_model$sample(
      data = stan_data,
      refresh = 200,
      chains = 2, iter_sampling = 500, iter_warmup = 500,
      parallel_chains = 2,
      adapt_delta = 0.8,
      max_treedepth = 10,
      show_exceptions = FALSE
    )

    summ <- stanfit$summary()
    fit_time <- round(stanfit$time()[["total"]] / 60, 1)
    cat("  Fit time:", fit_time, "minutes\n")

    stanfit$save_object(paste0(output_dir, "/", species_f, "_iCAR_New_",
                               firstYear, "_", lastYear, "_stanfit.rds"))
    saveRDS(summ, out_summ_file)

    ## --- Convergence check ---
    max_rhat <- max(summ$rhat, na.rm = TRUE)
    min_ess  <- min(summ$ess_bulk, na.rm = TRUE)
    cat("  Max Rhat:", round(max_rhat, 4),
        " | Min ESS:", round(min_ess, 0), "\n")

    ## --- Leave-future-out CV (final year only) ---
    ynext <- lastYear
    full_data <- new_data

    obs_df_fit <- full_data %>%
      filter(r_year <= ynext - 1) %>%
      mutate(observer = as.integer(factor(ObsN)))

    stan_data_cv <- list(
      count      = obs_df_fit$count,
      year       = obs_df_fit$year,
      route      = obs_df_fit$routeF,
      firstyr    = obs_df_fit$firstyr,
      observer   = obs_df_fit$observer,
      nobservers = max(obs_df_fit$observer),
      nyears     = max(obs_df_fit$year),
      nroutes    = nrow(route_map),
      ncounts    = length(obs_df_fit$count),
      fixedyear  = floor(max(obs_df_fit$year) / 2),
      N_edges    = car_stan_dat$N_edges,
      node1      = car_stan_dat$node1,
      node2      = car_stan_dat$node2,
      sd_alpha_prior = sd_alpha_prior
    )

    obs_df <- obs_df_fit %>% select(observer, ObsN) %>% distinct()

    obs_df_predict <- full_data %>%
      filter(r_year == ynext) %>%
      left_join(., obs_df, by = "ObsN") %>%
      mutate(observer = ifelse(!is.na(observer), observer, 0))

    stan_data_cv[["route_pred"]]    <- obs_df_predict$routeF
    stan_data_cv[["count_pred"]]    <- obs_df_predict$count
    stan_data_cv[["firstyr_pred"]]  <- obs_df_predict$firstyr
    stan_data_cv[["observer_pred"]] <- obs_df_predict$observer
    stan_data_cv[["ncounts_pred"]]  <- length(obs_df_predict$count)

    cv_lppd <- NA
    cv_cor  <- NA

    if (stan_data_cv[["ncounts_pred"]] > 0) {
      cv_fit <- cv_model$sample(
        data = stan_data_cv,
        refresh = 100,
        chains = 2, iter_sampling = 500, iter_warmup = 500,
        parallel_chains = 2,
        adapt_delta = 0.8,
        max_treedepth = 10,
        show_exceptions = FALSE
      )

      log_lik_full <- posterior_samples(fit = cv_fit, parm = "log_lik", dims = "i")
      log_lik_summ <- posterior_sums(log_lik_full, quantiles = NULL, dims = "i")
      names(log_lik_summ) <- paste0("log_lik_", names(log_lik_summ))

      E_pred_full <- posterior_samples(fit = cv_fit, parm = "E_pred", dims = "i")
      E_pred_summ <- posterior_sums(E_pred_full, quantiles = NULL, dims = "i")
      names(E_pred_summ) <- paste0("E_pred_", names(E_pred_summ))

      cv_results <- bind_cols(obs_df_predict, log_lik_summ, E_pred_summ)
      cv_results$E_pred_count <- exp(cv_results$E_pred_mean)

      cv_lppd <- mean(cv_results$log_lik_mean)
      cv_cor  <- cor(cv_results$count, cv_results$E_pred_count, method = "spearman")

      cat("  CV lppd:", round(cv_lppd, 4),
          " | Spearman r:", round(cv_cor, 3), "\n")
    } else {
      cat("  No prediction data for year", ynext, "— CV skipped.\n")
    }

    ## --- Extract trends and build trend map ---
    beta_samples <- posterior_samples(fit = stanfit, parm = "beta", dims = "route")
    beta_summ <- posterior_sums(beta_samples, dims = "route") %>%
      mutate(trend     = 100 * (exp(mean) - 1),
             trend_lci = 100 * (exp(lci) - 1),
             trend_uci = 100 * (exp(uci) - 1),
             routeF    = route)

    alpha_samples <- posterior_samples(fit = stanfit, parm = "alpha", dims = "route")
    alpha_summ <- posterior_sums(alpha_samples, dims = "route") %>%
      mutate(rel_abundance = exp(mean), routeF = route)

    trend_map_data <- route_map %>%
      inner_join(beta_summ, by = "routeF") %>%
      inner_join(alpha_summ %>% select(routeF, rel_abundance), by = "routeF")

    trend_breaks  <- c(-Inf, -4, -2, -0.5, 0.5, 2, 4, Inf)
    trend_labels  <- c("< -4%", "-4 to -2%", "-2 to -0.5%", "-0.5 to 0.5%",
                        "0.5 to 2%", "2 to 4%", "> 4%")
    trend_colours <- c("#d73027", "#fc8d59", "#fee08b", "#ffffbf",
                       "#d9ef8b", "#91cf60", "#1a9850")

    trend_map_data$trend_cat <- cut(trend_map_data$trend,
                                    breaks = trend_breaks, labels = trend_labels)

    # Trend map
    p_trend <- ggplot() +
      geom_sf(data = realized_strata_map, fill = grey(0.95), colour = grey(0.8)) +
      geom_sf(data = trend_map_data,
              aes(colour = trend_cat, size = rel_abundance), alpha = 0.8) +
      scale_colour_manual(values = trend_colours, name = "Annual trend", drop = FALSE) +
      scale_size_continuous(range = c(0.5, 4), name = "Relative\nabundance") +
      labs(title = paste(species, "— iCAR route-level trends", firstYear, "-", lastYear),
           subtitle = paste("Routes:", max(new_data$routeF),
                            " | Median trend:", round(median(trend_map_data$trend), 2), "%")) +
      theme_bw() +
      theme(legend.position = "right",
            text = element_text(family = "serif", size = 11))

    pdf(paste0("figures/", species_f, "_iCAR_trend_map_",
               firstYear, "_", lastYear, ".pdf"),
        width = 11, height = 8.5)
    print(p_trend)
    dev.off()

    png(paste0("figures/", species_f, "_iCAR_trend_map_",
               firstYear, "_", lastYear, ".png"),
        width = 11, height = 8.5, units = "in", res = 300)
    print(p_trend)
    dev.off()

    ## --- Collect summary row ---
    sp_result <- data.frame(
      species     = species,
      code        = grassland_spp$Code[grassland_spp$Common.Name == species][1],
      n_routes    = max(new_data$routeF),
      n_obs       = nrow(new_data),
      median_trend = round(median(trend_map_data$trend), 3),
      n_declining = sum(trend_map_data$trend < -0.5),
      n_increasing = sum(trend_map_data$trend > 0.5),
      n_stable    = sum(trend_map_data$trend >= -0.5 & trend_map_data$trend <= 0.5),
      max_rhat    = round(max_rhat, 4),
      min_ess     = round(min_ess, 0),
      cv_lppd     = round(cv_lppd, 4),
      cv_spearman = round(cv_cor, 3),
      fit_minutes = fit_time,
      stringsAsFactors = FALSE
    )

    all_results <- bind_rows(all_results, sp_result)

    cat("  Median trend:", round(median(trend_map_data$trend), 2), "% per year\n")
    cat("  Done.\n")

  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    cat("  Skipping", species, "\n")
  })

}  # end species loop

## =============================================================================
## Summary table across all species
## =============================================================================

cat("\n\n================================================================\n")
cat("  ALL SPECIES SUMMARY\n")
cat("================================================================\n\n")

if (!is.null(all_results) && nrow(all_results) > 0) {
  print(as.data.frame(all_results), row.names = FALSE)

  saveRDS(all_results, paste0(output_dir, "/grassland_iCAR_summary_",
                               firstYear, "_", lastYear, ".rds"))
  write.csv(all_results, paste0(output_dir, "/grassland_iCAR_summary_",
                                 firstYear, "_", lastYear, ".csv"),
            row.names = FALSE)

  # Summary visualization
  p_summary <- ggplot(all_results,
                       aes(x = reorder(species, median_trend),
                           y = median_trend)) +
    geom_col(aes(fill = median_trend > 0), alpha = 0.8, show.legend = FALSE) +
    scale_fill_manual(values = c("TRUE" = "#1a9850", "FALSE" = "#d73027")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    labs(title = paste("Grassland bird trends — iCAR model", firstYear, "-", lastYear),
         subtitle = paste(nrow(all_results), "species"),
         x = "", y = "Median annual % trend") +
    theme_bw() +
    theme(text = element_text(family = "serif", size = 12))

  pdf(paste0("figures/grassland_species_trends_summary_",
             firstYear, "_", lastYear, ".pdf"),
      width = 10, height = max(6, nrow(all_results) * 0.35))
  print(p_summary)
  dev.off()

  png(paste0("figures/grassland_species_trends_summary_",
             firstYear, "_", lastYear, ".png"),
      width = 10, height = max(6, nrow(all_results) * 0.35),
      units = "in", res = 300)
  print(p_summary)
  dev.off()

  cat("\nSummary table and figure saved.\n")
} else {
  cat("No species were successfully fitted.\n")
}

cat("\nDone.\n")
