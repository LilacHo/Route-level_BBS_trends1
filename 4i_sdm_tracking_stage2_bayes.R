# 4i_sdm_tracking_stage2_bayes.R
#
# Stage 2 of the Bayesian "SDM tracking" analysis: a hierarchical Bayesian
# random-effects meta-analysis of the per-species quantities produced by stage 1
# (4h_sdm_tracking_stage1_draws.R). Fast: one row per species. Runs anywhere.
#
# Question: does SDM-projected habitat-suitability change track observed
# population trends, and for which species / groups does it not?
#
# Main quantity: each species' TRACKING SLOPE (%/yr per SDM class step;
# positive = trends better where suitability is projected to improve).
# Each species contributes its stage-1 posterior mean (y) and SD (se); the
# model (models/sdm_tracking_meta.stan) is
#   y[s] ~ Normal(theta[s], se[s]),  theta[s] ~ Normal(mu + group effects, tau)
# with group effects partially pooled. It reports:
#   mu    average species' tracking slope (P(mu > 0) = probability of tracking)
#   tau   between-species SD of true slopes, net of estimation error
#   sigma for each grouping: how much groups truly differ (near 0 = no difference)
#   theta each species' shrunken slope -> classified as
#         Tracks       P(theta > 0) >= class_threshold
#         Contrary     P(theta < 0) >= class_threshold
#         Inconclusive otherwise
#
# Models fitted (all on the same species, so LOO comparisons are valid):
#   tracking slope:  overall | overall (Student-t) | overall (Student-t, species-level; used to classify species) |
#                    habitat group | family |
#                    order + family | diet (8) | migration, Red List (3), invasive (US-RIIS): fixed effects
#   worsening-side slope (classes 1-4) and improving-side slope (classes 4-7): overall
#   Contraction - Expansion: overall  (for comparison with the glmmTMB analysis, 4f)
#
# Reads:  output/files/sdm_tracking/sdm_tracking_stage1_species_<tag>_<scenario>_<firstYear>_<lastYear>.csv
#         data/spp_names_codes_group_aou_taxa_redlist_diet_migration_invasive.csv  (traits)
#         models/sdm_tracking_meta.stan
# Writes (output/files/sdm_tracking/):
#   sdm_tracking_stage2_SUMMARY_<suffix>.txt   <- the file to read / share; all key numbers
#   sdm_tracking_stage2_models_<suffix>.csv    one row per model: mu, tau, sigmas, LOO, diagnostics
#   sdm_tracking_stage2_groups_<suffix>.csv    one row per model x group level
#   sdm_tracking_stage2_species_<suffix>.csv   one row per species: stage-1 and shrunken slopes, class
#
# Needs: tidyverse, cmdstanr (with CmdStan installed), posterior, loo.

library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(loo)

here::i_am("4i_sdm_tracking_stage2_bayes.R")

# Settings -------------------------------------------------------------------
model_tag <- "anthro"; scenario <- "rcp45"; firstYear <- 2010; lastYear <- 2025
# Must match exclude_classes in 4h_sdm_tracking_stage1_draws.R: integer(0) = main analysis,
# c(7) = colonization excluded (reads and writes files tagged "_no7").
exclude_classes <- c(7)
run_tag <- if (length(exclude_classes) > 0) paste0("_no", paste(sort(exclude_classes), collapse = "")) else ""
suffix <- paste0(paste(model_tag, scenario, firstYear, lastYear, sep = "_"), run_tag)
stage1_csv <- here::here("output", "files", "sdm_tracking", paste0("sdm_tracking_stage1_species_", suffix, ".csv"))
traits_csv <- here::here("data", "spp_names_codes_group_aou_taxa_redlist_diet_migration_invasive.csv")
stan_file  <- here::here("models", "sdm_tracking_meta.stan")
out_dir    <- here::here("output", "files", "sdm_tracking")

class_threshold <- 0.90   # posterior probability needed to call a species/group Tracks or Contrary
prior_scale_slope <- 1    # %/yr per class step (observed scale ~0.2)
prior_scale_ce    <- 3    # %/yr (Contraction - Expansion, observed scale ~1)
chains <- 4; iter_warmup <- 1000; iter_sampling <- 2000; adapt_delta <- 0.99
seed <- 20261014

out_txt     <- file.path(out_dir, paste0("sdm_tracking_stage2_SUMMARY_", suffix, ".txt"))
out_models  <- file.path(out_dir, paste0("sdm_tracking_stage2_models_", suffix, ".csv"))
out_groups  <- file.path(out_dir, paste0("sdm_tracking_stage2_groups_", suffix, ".csv"))
out_species <- file.path(out_dir, paste0("sdm_tracking_stage2_species_", suffix, ".csv"))

# Data -----------------------------------------------------------------------
stage1 <- read.csv(stage1_csv, stringsAsFactors = FALSE)
traits <- read.csv(traits_csv, stringsAsFactors = FALSE) %>%
  filter(in_bbs %in% TRUE) %>% distinct(Code, .keep_all = TRUE) %>%
  transmute(species_code = Code, group = Group, order, family,
            diet = case_when(diet %in% c("Seed", "Plant", "Herbivore") ~ "Plant-eater",
                             diet %in% c("Vertebrate", "Carnivore", "Scavenger") ~ "Meat-eater", TRUE ~ diet),
            migration,
            redlist = case_when(redlistCategory %in% c("Vulnerable", "Endangered", "Critically Endangered") ~ "Threatened (VU + EN)",
                                TRUE ~ redlistCategory),
            invasive = if_else(invasive == "Not listed", "Native (not listed)", invasive))
d_all <- stage1 %>% left_join(traits, by = "species_code")

mod <- cmdstan_model(stan_file)

# Helpers --------------------------------------------------------------------
idx <- function(v) { f <- factor(v); list(i = as.integer(f), K = nlevels(f), levels = levels(f)) }

# explicit = 1: species-level Student-t hierarchy (for classification; not used in LOO)
# fixed1   = 1: grouping 1 as separate effects without a between-group SD (2-3 level traits)
fit_meta <- function(d, ycol, secol, grp1 = NULL, grp2 = NULL, student = 0, prior_scale,
                     explicit = 0, fixed1 = 0) {
  g1 <- if (is.null(grp1)) list(i = rep(0L, nrow(d)), K = 0L, levels = character()) else idx(d[[grp1]])
  g2 <- if (is.null(grp2)) list(i = rep(0L, nrow(d)), K = 0L, levels = character()) else idx(d[[grp2]])
  run <- function() mod$sample(data = list(N = nrow(d), y = d[[ycol]], se = d[[secol]],
                                           K1 = g1$K, g1 = g1$i, K2 = g2$K, g2 = g2$i,
                                           student = student, explicit_t = explicit, fixed1 = fixed1,
                                           prior_scale = prior_scale),
                               chains = chains, parallel_chains = chains, iter_warmup = iter_warmup,
                               iter_sampling = iter_sampling, adapt_delta = adapt_delta, seed = seed,
                               refresh = 0, show_messages = FALSE)
  fit <- run()
  # On Windows a freshly compiled model executable can be blocked on its first launch
  # (e.g. by antivirus scanning), so that no chain starts. Retry once after a pause.
  if (all(fit$return_codes() != 0)) {
    message("No chains finished; retrying once")
    Sys.sleep(10)
    fit <- run()
    if (all(fit$return_codes() != 0)) stop("No chains finished on retry; rerun with show_messages = TRUE to see why")
  }
  list(fit = fit, g1 = g1, g2 = g2, d = d, ycol = ycol, explicit = explicit, fixed1 = fixed1)
}

q <- function(v) c(mean = mean(v), lo = unname(quantile(v, 0.025)), hi = unname(quantile(v, 0.975)))
fmt_ci <- function(v, digits = 2) { z <- q(v); sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), z["mean"], z["lo"], z["hi"]) }
classify <- function(p_pos) case_when(p_pos >= class_threshold ~ "Tracks",
                                      1 - p_pos >= class_threshold ~ "Contrary",
                                      TRUE ~ "Inconclusive")

summarise_model <- function(m, name, response) {
  fit <- m$fit
  dr <- fit$draws(c("mu", "tau", "sigma1", "sigma2", "nu"), format = "draws_df")
  pooled1 <- m$g1$K > 0 && m$fixed1 == 0
  sum_vars <- c("mu", "tau", if (m$g1$K > 0) "z1", if (pooled1) "sigma1", if (m$g2$K > 0) c("sigma2", "z2"),
                if (all(is.finite(dr$nu))) "nu")
  sm <- fit$summary(sum_vars) %>% filter(!is.na(rhat))
  diag <- fit$diagnostic_summary(quiet = TRUE)
  lo <- suppressWarnings(fit$loo(variables = "log_lik", cores = 1))
  # share of species whose true value is negative: Normal or Student-t species distribution
  neg <- if (all(is.finite(dr$nu))) pt(-dr$mu / dr$tau, df = dr$nu) else pnorm(-dr$mu / dr$tau)
  row <- tibble(model = name, response = response, n_species = nrow(m$d),
                mu_mean = mean(dr$mu), mu_lo = q(dr$mu)[["lo"]], mu_hi = q(dr$mu)[["hi"]],
                p_mu_pos = mean(dr$mu > 0),
                tau_mean = mean(dr$tau), tau_lo = q(dr$tau)[["lo"]], tau_hi = q(dr$tau)[["hi"]],
                pct_species_negative = 100 * mean(neg),
                group_effects = if (m$g1$K == 0) "none" else if (m$fixed1 == 1) "fixed" else "partially pooled",
                sigma1_mean = if (pooled1) mean(dr$sigma1) else NA, sigma1_lo = if (pooled1) q(dr$sigma1)[["lo"]] else NA,
                sigma1_hi = if (pooled1) q(dr$sigma1)[["hi"]] else NA,
                sigma2_mean = if (m$g2$K > 0) mean(dr$sigma2) else NA, sigma2_lo = if (m$g2$K > 0) q(dr$sigma2)[["lo"]] else NA,
                sigma2_hi = if (m$g2$K > 0) q(dr$sigma2)[["hi"]] else NA,
                nu_mean = if (all(is.finite(dr$nu))) mean(dr$nu) else NA,
                elpd_loo = lo$estimates["elpd_loo", "Estimate"], elpd_loo_se = lo$estimates["elpd_loo", "SE"],
                pareto_k_over_0.7 = sum(lo$diagnostics$pareto_k > 0.7),
                max_rhat = max(sm$rhat, na.rm = TRUE), min_ess_bulk = min(sm$ess_bulk, na.rm = TRUE),
                divergences = sum(diag$num_divergent))
  # posterior predictive check: does the spread of simulated species values match the observed?
  yrep <- fit$draws("y_rep", format = "draws_matrix")
  sd_rep <- apply(yrep, 1, sd)
  row$ppc_p_sd <- mean(sd_rep >= sd(m$d[[m$ycol]]))
  list(row = row, loo = lo, draws = dr)
}

group_table <- function(m, name, which = 1) {
  g <- if (which == 1) m$g1 else m$g2
  if (g$K == 0) return(NULL)
  dr <- m$fit$draws(c("mu", paste0("a", which)), format = "draws_matrix")
  mu <- dr[, "mu"]
  map_dfr(seq_len(g$K), function(k) {
    v <- as.vector(mu + dr[, paste0("a", which, "[", k, "]")])
    members <- m$d[g$i == k, ]
    tibble(model = name, grouping = if (which == 1) "g1" else "g2", level = g$levels[k],
           n_species = nrow(members),
           mean = mean(v), lo = q(v)[["lo"]], hi = q(v)[["hi"]], p_pos = mean(v > 0),
           class = classify(mean(v > 0)))
  })
}

species_theta <- function(m) {
  th <- m$fit$draws("theta", format = "draws_matrix")
  tibble(species_code = m$d$species_code,
         theta_mean = colMeans(th),
         theta_lo = apply(th, 2, quantile, 0.025), theta_hi = apply(th, 2, quantile, 0.975),
         theta_p_pos = colMeans(th > 0))
}

# Fit ------------------------------------------------------------------------
d <- d_all %>% filter(status == "ok", is.finite(slope_mean), slope_sd > 0)
cat("Species with a tracking slope:", nrow(d), "of", nrow(d_all), "\n")
# name, grouping 1, grouping 2, student, explicit species-level t, fixed group effects
t_name <- "Overall (Student-t, species-level)"
specs <- list(
  list("Overall",             NULL,        NULL,     0, 0, 0),
  list("Overall (Student-t)", NULL,        NULL,     1, 0, 0),
  list(t_name,                NULL,        NULL,     1, 1, 0),
  list("Habitat group",       "group",     NULL,     0, 0, 0),
  list("Family",              "family",    NULL,     0, 0, 0),
  list("Order + family",      "order",     "family", 0, 0, 0),
  list("Diet (8)",            "diet",      NULL,     0, 0, 0),
  list("Migration",           "migration", NULL,     0, 0, 1),
  list("Red List (3)",        "redlist",   NULL,     0, 0, 1),
  list("Invasive (US-RIIS)",  "invasive",  NULL,     0, 0, 1)
)
fits <- list(); rows <- list(); loos <- list(); groups <- list()
for (s in specs) {
  cat("Fitting:", s[[1]], "\n")
  m <- fit_meta(d, "slope_mean", "slope_sd", s[[2]], s[[3]], s[[4]], prior_scale_slope,
                explicit = s[[5]], fixed1 = s[[6]])
  r <- summarise_model(m, s[[1]], "Tracking slope")
  if (s[[5]] == 1) {   # conditional log_lik: not comparable in LOO
    r$row$elpd_loo <- NA; r$row$elpd_loo_se <- NA; r$row$pareto_k_over_0.7 <- NA
  } else {
    loos[[s[[1]]]] <- r$loo
  }
  fits[[s[[1]]]] <- m; rows[[s[[1]]]] <- r$row
  groups[[s[[1]]]] <- bind_rows(group_table(m, s[[1]], 1), group_table(m, s[[1]], 2))
}
side <- list()
for (resp in c("slope_low", "slope_up", "ce")) {
  dd <- d_all %>% filter(status == "ok", is.finite(.data[[paste0(resp, "_mean")]]), .data[[paste0(resp, "_sd")]] > 0)
  top_class <- max(setdiff(1:7, exclude_classes))
  lab <- c(slope_low = "Worsening-side slope (classes 1-4)",
           slope_up = paste0("Improving-side slope (classes 4-", top_class, ")"),
           ce = "Contraction - Expansion")[[resp]]
  cat("Fitting:", lab, "\n")
  m <- fit_meta(dd, paste0(resp, "_mean"), paste0(resp, "_sd"), NULL, NULL, 0,
                if (resp == "ce") prior_scale_ce else prior_scale_slope)
  side[[resp]] <- summarise_model(m, "Overall", lab)$row
}

# unname(): a named list of one-row tibbles would otherwise be read by bind_rows() as one data frame
side_tab   <- bind_rows(unname(side))
models_out <- bind_rows(c(unname(rows), unname(side)))
# Each model vs the no-trait Normal "Overall" model: positive elpd_diff = adding the
# trait (or Student-t errors) improves leave-one-species-out prediction.
pair_diff <- function(nm) {
  if (nm == "Overall") return(c(0, 0))
  d_i <- loos[[nm]]$pointwise[, "elpd_loo"] - loos[["Overall"]]$pointwise[, "elpd_loo"]
  c(sum(d_i), sqrt(length(d_i)) * sd(d_i))
}
cmp_tab <- map_dfr(names(loos), function(nm) {
  v <- pair_diff(nm)
  tibble(model = nm, response = "Tracking slope", elpd_diff = v[1], se_diff = v[2])
})
models_out <- models_out %>% left_join(cmp_tab, by = c("model", "response"))
groups_out <- bind_rows(unname(groups))

sp0 <- species_theta(fits[["Overall"]]) %>% rename_with(~ paste0(.x, "_overall"), -species_code)
spt <- species_theta(fits[[t_name]]) %>% rename_with(~ paste0(.x, "_studentt"), -species_code)
spg <- species_theta(fits[["Habitat group"]]) %>% rename_with(~ paste0(.x, "_groupmodel"), -species_code)
species_out <- d %>%
  select(species, species_code, group, family, order, diet, migration, redlist, invasive, n_routes, n_classes,
         starts_with("n_class"), slope_mean, slope_sd, slope_p_pos, slope_low_mean, slope_up_mean, ce_mean) %>%
  left_join(sp0, by = "species_code") %>% left_join(spt, by = "species_code") %>% left_join(spg, by = "species_code") %>%
  mutate(class_stage1 = classify(slope_p_pos),
         class_overall = classify(theta_p_pos_overall),
         class_studentt = classify(theta_p_pos_studentt),
         class_groupmodel = classify(theta_p_pos_groupmodel)) %>%
  arrange(theta_mean_studentt)

write.csv(models_out, out_models, row.names = FALSE)
write.csv(groups_out, out_groups, row.names = FALSE)
write.csv(species_out, out_species, row.names = FALSE)

# Plain-text summary ---------------------------------------------------------
con <- file(out_txt, open = "wt", encoding = "UTF-8")
w <- function(...) { cat(..., "\n", sep = "", file = con); cat(..., "\n", sep = "") }
f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
mr <- function(name) models_out %>% filter(model == name, response == "Tracking slope")

w("SDM TRACKING - STAGE 2 SUMMARY (", suffix, ")   ", format(Sys.time(), "%Y-%m-%d %H:%M"))
w("Tracking slope = change in route trend (%/yr) per SDM class step; positive = tracks the projection.")
w("Species with a slope: ", nrow(d), " (stage-1 rows: ", nrow(d_all), "; status counts: ",
  paste(names(table(d_all$status)), table(d_all$status), sep = " = ", collapse = "; "), ")")
w("Class threshold: P >= ", class_threshold, ".  Intervals are 95% credible intervals.")
w("")
w("== 1. DOES SDM CHANGE TRACK OBSERVED TRENDS? (overall model) ==")
o <- mr("Overall")
w("mu (average species' slope): ", f2(o$mu_mean, 3), " [", f2(o$mu_lo, 3), ", ", f2(o$mu_hi, 3), "]   P(mu > 0) = ", f2(o$p_mu_pos, 4))
w("tau (between-species SD of true slopes): ", f2(o$tau_mean, 3), " [", f2(o$tau_lo, 3), ", ", f2(o$tau_hi, 3), "]")
w("Expected % of species with a negative true slope: ", f2(o$pct_species_negative, 1), "%")
ot <- mr("Overall (Student-t)")
w("Student-t version (marginal, used for LOO): mu ", f2(ot$mu_mean, 3), " [", f2(ot$mu_lo, 3), ", ", f2(ot$mu_hi, 3), "], nu = ", f2(ot$nu_mean, 1),
  ";  elpd vs Normal: ", f2(ot$elpd_loo - o$elpd_loo, 1))
os <- mr(t_name)
w("Student-t species-level (used for species classification): mu ", f2(os$mu_mean, 3), " [", f2(os$mu_lo, 3), ", ", f2(os$mu_hi, 3),
  "], tau (scale) ", f2(os$tau_mean, 3), ", nu = ", f2(os$nu_mean, 1), ", expected % species negative ", f2(os$pct_species_negative, 1), "%",
  ";  divergences ", os$divergences)
w("Raw stage-1 (unpooled) slopes: median ", f2(median(d$slope_mean), 3), "; ", sum(d$slope_p_pos >= class_threshold), " species P(>0) >= ",
  class_threshold, ", ", sum(1 - d$slope_p_pos >= class_threshold), " species P(<0) >= ", class_threshold)
w("")
w("== 2. WORSENING vs IMPROVING SIDE, and CONTRACTION - EXPANSION ==")
for (r in seq_len(nrow(side_tab))) {
  x <- side_tab[r, ]
  w(sprintf("%-38s n = %3d  mu %s  P(mu>0) = %.3f  tau %.3f", x$response, x$n_species,
            paste0(f2(x$mu_mean, 3), " [", f2(x$mu_lo, 3), ", ", f2(x$mu_hi, 3), "]"), x$p_mu_pos, x$tau_mean))
}
w("")
w("== 3. SPECIES CLASSIFICATION (shrunken slopes) ==")
for (cl in c("class_stage1", "class_overall", "class_studentt", "class_groupmodel")) {
  tb <- table(factor(species_out[[cl]], levels = c("Tracks", "Inconclusive", "Contrary")))
  w(sprintf("%-18s Tracks %3d | Inconclusive %3d | Contrary %3d", sub("class_", "", cl), tb[1], tb[2], tb[3]))
}
w("Contrary species (Student-t species-level model; Normal-model class in brackets): ")
ctr <- species_out %>% filter(class_studentt == "Contrary") %>% arrange(desc(theta_mean_studentt))
if (nrow(ctr) == 0) w("  none") else for (i in seq_len(nrow(ctr))) {
  x <- ctr[i, ]
  w(sprintf("  %-30s %-20s slope %s  P(>0)=%.3f  routes %4d  [%s]", x$species, x$group,
            paste0(f2(x$theta_mean_studentt, 3), " [", f2(x$theta_lo_studentt, 3), ", ", f2(x$theta_hi_studentt, 3), "]"),
            x$theta_p_pos_studentt, x$n_routes, x$class_overall))
}
w("Contrary under the Normal model only (not under Student-t): ",
  paste(species_out$species[species_out$class_overall == "Contrary" & species_out$class_studentt != "Contrary"], collapse = "; "))
w("Strongest-tracking 10 species (Student-t species-level model): ")
top <- species_out %>% arrange(desc(theta_mean_studentt)) %>% head(10)
for (i in seq_len(nrow(top))) w(sprintf("  %-30s %-20s slope %.3f  P(>0)=%.3f", top$species[i], top$group[i], top$theta_mean_studentt[i], top$theta_p_pos_studentt[i]))
w("")
w("== 4. WHICH TRAITS ORGANISE TRACKING? (model comparison) ==")
w("sigma = SD of true group effects (near 0 = groups do not differ). elpd_diff: LOO improvement over the no-trait Overall model (positive = trait helps; |diff| < ~2*se_diff = no clear difference).")
mt <- models_out %>% filter(response == "Tracking slope")
for (r in seq_len(nrow(mt))) {
  x <- mt[r, ]
  w(sprintf("%-36s sigma1 %-22s sigma2 %-22s elpd_diff %s  tau %.3f",
            x$model,
            if (x$group_effects == "fixed") "(fixed effects)" else if (is.na(x$sigma1_mean)) "-" else sprintf("%.3f [%.3f,%.3f]", x$sigma1_mean, x$sigma1_lo, x$sigma1_hi),
            if (is.na(x$sigma2_mean)) "-" else sprintf("%.3f [%.3f,%.3f]", x$sigma2_mean, x$sigma2_lo, x$sigma2_hi),
            if (is.na(x$elpd_diff)) "  n/a (not LOO-comparable)" else sprintf("%6.1f (se %4.1f)", x$elpd_diff, x$se_diff),
            x$tau_mean))
}
w("")
w("== 5. GROUP-LEVEL SLOPES (mu + group effect) ==")
for (nm in c("Habitat group", "Diet (8)", "Migration", "Red List (3)", "Invasive (US-RIIS)", "Family")) {
  gt <- groups_out %>% filter(model == nm, grouping == "g1") %>% arrange(desc(mean))
  if (nm == "Family") gt <- gt %>% filter(n_species >= 5)
  w("-- ", nm, if (nm == "Family") " (families with >= 5 species)" else "", " --")
  for (i in seq_len(nrow(gt))) {
    sc <- species_out %>% filter(.data[[c(`Habitat group` = "group", `Diet (8)` = "diet", Migration = "migration",
                                          `Red List (3)` = "redlist", `Invasive (US-RIIS)` = "invasive", Family = "family")[[nm]]]] == gt$level[i])
    tb <- table(factor(sc$class_studentt, levels = c("Tracks", "Inconclusive", "Contrary")))
    w(sprintf("  %-28s n=%3d  %s  P(>0)=%.3f  %-12s species T/I/C %d/%d/%d", gt$level[i], gt$n_species[i],
              paste0(f2(gt$mean[i], 3), " [", f2(gt$lo[i], 3), ", ", f2(gt$hi[i], 3), "]"), gt$p_pos[i], gt$class[i], tb[1], tb[2], tb[3]))
  }
}
w("")
w("== 6. DIAGNOSTICS ==")
for (r in seq_len(nrow(models_out))) {
  x <- models_out[r, ]
  w(sprintf("%-36s %-38s max Rhat %.3f  min ESS %5.0f  divergences %d  Pareto k>0.7: %d  PPC p(sd) %.2f",
            x$model, x$response, x$max_rhat, x$min_ess_bulk, x$divergences, x$pareto_k_over_0.7, x$ppc_p_sd))
}
close(con)
cat("\nWrote:\n ", out_txt, "\n ", out_models, "\n ", out_groups, "\n ", out_species, "\n")
