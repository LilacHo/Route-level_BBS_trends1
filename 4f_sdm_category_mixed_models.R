# 4f_sdm_category_mixed_models.R
#
# Tests whether SDM range-change category (Contraction / Stable / Expansion)
# predicts route-level population trend WITHIN species, with species as the
# unit of replication. This is the main analysis behind the manuscript's
# Results; robustness checks live in 4g_sdm_category_robustness.R.
#
# Why a mixed model: pooling all routes (4c's Kruskal-Wallis, 4d's category
# means) treats routes of the same species as independent and mixes the
# within-species effect with species composition (species doing well
# everywhere may simply have more Expansion routes). Here every species gets
# its own baseline trend (random intercept) and its own Contraction /
# Expansion effects (random slopes), so the fixed effects are the
# average-species within-species effects.
#
# Why Student-t errors: route trends are heavy-tailed (residual excess
# kurtosis ~6 under a Gaussian fit), so a normal residual lets a few extreme
# routes pull the estimates. A Student-t residual (df estimated) down-weights
# them. The tail check below confirms the fitted t captures the tails.
#
# Headline contrast: CONTRACTION - EXPANSION (C - E, %/yr), predicted to be
# NEGATIVE (routes in a species' Contraction areas trend worse than its
# Expansion routes).
#
# 1. Main model (glmmTMB, maximum likelihood):
#      trend ~ cat3 + (1 + cat3 | species_code), family = t_family()
#    cat3 reference level = Stable. C - E = beta_C - beta_E, SE from the
#    fixed-effect covariance (L = (0, 1, -1)); Wald z p-values.
#    Also: species random-effect SDs / correlations, the expected share of
#    species with a reversed (positive) C - E effect, species-level fitted
#    values (fixed + random), a residual tail check against the fitted t,
#    and within-species R2 of category (model-free).
# 2. Bird group (arctic dropped: 1 species, 2 Stable routes, no contrast):
#      trend ~ 0 + group + group:cat3 + (1 + cat3 | species_code), t errors
#    Per-group C - E with Benjamini-Hochberg q across groups, and a
#    likelihood-ratio test against group-specific baselines with shared
#    category effects: trend ~ 0 + group + cat3 + (1 + cat3 | species_code).
#
# Reads:  output/species_routes_covariates/per_species_sdm/*_route_trends_sdm.csv
#         (written by 3c_add_SDM_covariates.R)
# Writes: output/files/sdm_category_mixed_model_summary_2010_2025.csv
#           One row: data counts, fixed effects and C - E (SE, 95% CI, p),
#           t df / scale, random-effect SDs / correlations, expected % of
#           species reversed, tail check, within-species R2, group LRT and
#           species-level SDs with shared vs group-specific effects.
#         output/files/sdm_category_mixed_model_species_effects_2010_2025.csv
#           One row per species: route counts and raw means per category,
#           fitted Contraction / Stable / Expansion trends and C - E.
#         output/files/sdm_category_mixed_model_by_group_2010_2025.csv
#           One row per group: counts, Stable baseline, C - S, E - S and
#           C - E (SE, 95% CI, p, BH q).
# Write-ups (hand-built HTML, not generated here):
#   output/reports/sdm_category_trend_report_student_t_anthro_rcp45_2010_2025.html

library(here)
library(tidyverse)
library(glmmTMB)

here::i_am("4f_sdm_category_mixed_models.R")

# Settings -------------------------------------------------------------------
model_tag  <- "anthro"   # trend model: "anthro" or "base"
scenario   <- "rcp45"    # SDM scenario column: "rcp45" or "rcp85"
first_year <- 2010
last_year  <- 2025
require_route_converged <- TRUE   # same rule as 4c
excluded_groups <- c("arctic")    # no Contraction/Expansion routes at all

in_dir  <- here::here("output", "species_routes_covariates", "per_species_sdm")
out_dir <- here::here("output", "files")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
suffix <- paste0("_", first_year, "_", last_year, ".csv")
out_summary <- file.path(out_dir, paste0("sdm_category_mixed_model_summary", suffix))
out_species <- file.path(out_dir, paste0("sdm_category_mixed_model_species_effects", suffix))
out_group   <- file.path(out_dir, paste0("sdm_category_mixed_model_by_group", suffix))

z975 <- qnorm(0.975)
p_norm <- function(z) 2 * pnorm(-abs(z))

# Linear contrast of fixed effects: estimate and SE.
contrast <- function(b, V, l) {
  l <- l[names(b)]
  l[is.na(l)] <- 0
  c(est = sum(l * b), se = sqrt(drop(t(l) %*% V %*% l)))
}
# Weights for Contraction - Expansion.
ce_weights <- function(con_name, exp_name) setNames(c(1, -1), c(con_name, exp_name))
# glmmTMB's fixef()/vcov() return lists; take the conditional model.
tmb_b <- function(f) fixef(f)$cond
tmb_V <- function(f) as.matrix(vcov(f)$cond)
tmb_t_df <- function(f) {
  tryCatch(unname(family_params(f)[1]), error = function(e) exp(unname(f$fit$par["psi"])))
}
check_fit <- function(f, label) {
  ok <- f$fit$convergence == 0 && isTRUE(f$sdr$pdHess)
  if (!ok) warning(label, ": optimiser did not converge or Hessian not positive definite")
  ok
}

# Read ----------------------------------------------------------------------
sdm_files <- list.files(in_dir, pattern = "_route_trends_sdm\\.csv$", full.names = TRUE)
if (length(sdm_files) == 0) stop("No *_route_trends_sdm.csv files in ", in_dir, " -- run 3c first.")
all_routes <- sdm_files %>%
  map_dfr(~ read.csv(.x, stringsAsFactors = FALSE)) %>%
  as_tibble()
if (require_route_converged) {
  n_before <- nrow(all_routes)
  all_routes <- all_routes %>% filter(route_converged %in% TRUE)
  cat("require_route_converged = TRUE -> dropped", n_before - nrow(all_routes), "of", n_before, "rows\n")
}

x <- all_routes %>%
  filter(model == model_tag) %>%
  transmute(species, species_code, group, trend, cat7 = .data[[scenario]]) %>%
  filter(!is.na(cat7), cat7 %in% 1:7) %>%
  mutate(cat3 = factor(case_when(cat7 <= 3 ~ "Contraction", cat7 == 4 ~ "Stable", TRUE ~ "Expansion"),
                       levels = c("Stable", "Contraction", "Expansion")))
n_cat_by_sp <- x %>% group_by(species_code) %>% summarise(n_cat = n_distinct(cat3), .groups = "drop")
cat("\n==================", model_tag, scenario, "==================\n")
cat("routes:", nrow(x), " species:", n_distinct(x$species_code),
    " species with >=2 categories:", sum(n_cat_by_sp$n_cat >= 2), "\n")

# 1. Main model ----------------------------------------------------------------
fit <- glmmTMB(trend ~ cat3 + (1 + cat3 | species_code), data = x, family = t_family())
fit_ok <- check_fit(fit, "main model")

b <- tmb_b(fit); V <- tmb_V(fit)
fe <- coef(summary(fit))$cond
g <- contrast(b, V, ce_weights("cat3Contraction", "cat3Expansion"))
t_df <- tmb_t_df(fit)
cat(sprintf("C - E %.2f (SE %.2f, 95%% CI %.2f to %.2f, p = %.2g); t df = %.2f\n",
            g[["est"]], g[["se"]], g[["est"]] - z975 * g[["se"]], g[["est"]] + z975 * g[["se"]],
            p_norm(g[["est"]] / g[["se"]]), t_df))

# Between-species variation in the category effects
vc <- VarCorr(fit)$cond$species_code
sd_re <- attr(vc, "stddev")
cor_re <- attr(vc, "correlation")
sd_ce <- sqrt(vc["cat3Contraction", "cat3Contraction"] + vc["cat3Expansion", "cat3Expansion"] -
                2 * vc["cat3Contraction", "cat3Expansion"])

# Tail check: share of standardised residuals beyond the fitted t's 1% and
# 0.1% two-sided cutoffs -- should be close to 1% and 0.1% if the t fits.
r_t <- residuals(fit) / sigma(fit)
pct_beyond_1pct   <- 100 * mean(abs(r_t) > qt(0.995, t_df))
pct_beyond_0.1pct <- 100 * mean(abs(r_t) > qt(0.9995, t_df))

# Predictive strength (model-free): share of within-species variation in
# trend explained by category, vs variation explained by species identity.
x_dm <- x %>% left_join(n_cat_by_sp, by = "species_code") %>% filter(n_cat >= 2) %>%
  group_by(species_code) %>% mutate(trend_dm = trend - mean(trend)) %>% ungroup()
r2_within  <- summary(lm(trend_dm ~ cat3, data = x_dm))$r.squared
r2_species <- summary(lm(trend ~ species_code, data = x))$r.squared

# Species-level estimates (fixed + random effects)
sp_raw <- x %>%
  group_by(species_code) %>%
  summarise(species = first(species), group = first(group), n_routes = n(),
            n_con = sum(cat3 == "Contraction"), n_stable = sum(cat3 == "Stable"),
            n_exp = sum(cat3 == "Expansion"),
            raw_con = if (n_con > 0) mean(trend[cat3 == "Contraction"]) else NA_real_,
            raw_stable = if (n_stable > 0) mean(trend[cat3 == "Stable"]) else NA_real_,
            raw_exp = if (n_exp > 0) mean(trend[cat3 == "Expansion"]) else NA_real_,
            .groups = "drop")
sp_fit <- coef(fit)$cond$species_code %>%
  rownames_to_column("species_code") %>%
  as_tibble() %>%
  transmute(species_code,
            fit_stable = `(Intercept)`,
            fit_con = `(Intercept)` + cat3Contraction,
            fit_exp = `(Intercept)` + cat3Expansion,
            fit_con_minus_exp = cat3Contraction - cat3Expansion)
species_out <- sp_raw %>%
  left_join(sp_fit, by = "species_code") %>%
  mutate(model = model_tag, scenario = scenario, .before = 1)
both <- species_out %>% filter(n_con >= 1, n_exp >= 1)

# 2. Bird group ------------------------------------------------------------------
xg <- x %>% filter(!group %in% excluded_groups) %>% mutate(group = factor(group))
fit_g  <- glmmTMB(trend ~ 0 + group + group:cat3 + (1 + cat3 | species_code), data = xg, family = t_family())
fit_g0 <- glmmTMB(trend ~ 0 + group + cat3 + (1 + cat3 | species_code), data = xg, family = t_family())
invisible(check_fit(fit_g, "group model")); invisible(check_fit(fit_g0, "shared-effect group model"))
lrt <- anova(fit_g0, fit_g)   # glmmTMB: parameter count in "Df", test df in "Chi Df"

b_g <- tmb_b(fit_g); V_g <- tmb_V(fit_g)
group_tab <- map_dfr(levels(xg$group), function(gr) {
  nm_s <- paste0("group", gr)
  nm_c <- paste0("group", gr, ":cat3Contraction")
  nm_e <- paste0("group", gr, ":cat3Expansion")
  if (!all(c(nm_c, nm_e) %in% names(b_g))) return(NULL)  # group lacks a category
  ce <- contrast(b_g, V_g, ce_weights(nm_c, nm_e))
  tibble(group = gr, stable = b_g[[nm_s]],
         con_minus_stable = b_g[[nm_c]], con_minus_stable_se = sqrt(V_g[nm_c, nm_c]),
         exp_minus_stable = b_g[[nm_e]], exp_minus_stable_se = sqrt(V_g[nm_e, nm_e]),
         con_minus_exp = ce[["est"]], con_minus_exp_se = ce[["se"]])
}) %>%
  mutate(ci_lo = con_minus_exp - z975 * con_minus_exp_se,
         ci_hi = con_minus_exp + z975 * con_minus_exp_se,
         p = p_norm(con_minus_exp / con_minus_exp_se),
         q_bh = p.adjust(p, method = "BH"))
group_counts <- xg %>% group_by(group) %>%
  summarise(n_species = n_distinct(species_code), n_routes = n(),
            n_con = sum(cat3 == "Contraction"), n_stable = sum(cat3 == "Stable"),
            n_exp = sum(cat3 == "Expansion"), .groups = "drop") %>%
  mutate(group = as.character(group))
group_out <- group_tab %>%
  left_join(group_counts, by = "group") %>%
  mutate(model = model_tag, scenario = scenario, .before = 1) %>%
  arrange(con_minus_exp)

# How much of the between-species spread do groups explain?
sd_shared <- attr(VarCorr(fit_g0)$cond$species_code, "stddev")
sd_group  <- attr(VarCorr(fit_g)$cond$species_code, "stddev")

cat(sprintf("Groups differ: chisq = %.1f, df = %d, p = %.2g\n",
            lrt$Chisq[2], lrt$`Chi Df`[2], lrt$`Pr(>Chisq)`[2]))
print(as.data.frame(group_out %>% select(group, n_species, con_minus_exp, ci_lo, ci_hi, q_bh)),
      digits = 3, row.names = FALSE)

# Write ---------------------------------------------------------------------------
summary_out <- tibble(
  model = model_tag, scenario = scenario,
  n_routes = nrow(x), n_species = n_distinct(x$species_code),
  n_species_2cat = sum(n_cat_by_sp$n_cat >= 2),
  n_con = sum(x$cat3 == "Contraction"), n_stable = sum(x$cat3 == "Stable"), n_exp = sum(x$cat3 == "Expansion"),
  converged = fit_ok,
  stable_avg_species = fe["(Intercept)", 1], stable_avg_species_se = fe["(Intercept)", 2],
  con_minus_stable = fe["cat3Contraction", 1], con_minus_stable_se = fe["cat3Contraction", 2],
  con_minus_stable_p = fe["cat3Contraction", 4],
  exp_minus_stable = fe["cat3Expansion", 1], exp_minus_stable_se = fe["cat3Expansion", 2],
  exp_minus_stable_p = fe["cat3Expansion", 4],
  con_minus_exp = g[["est"]], con_minus_exp_se = g[["se"]],
  con_minus_exp_ci_lo = g[["est"]] - z975 * g[["se"]],
  con_minus_exp_ci_hi = g[["est"]] + z975 * g[["se"]],
  con_minus_exp_p = p_norm(g[["est"]] / g[["se"]]),
  t_df = t_df, t_sigma = sigma(fit),
  sd_species_intercept = sd_re[["(Intercept)"]], sd_species_con = sd_re[["cat3Contraction"]],
  sd_species_exp = sd_re[["cat3Expansion"]], sd_species_con_minus_exp = sd_ce,
  cor_intercept_con = cor_re["(Intercept)", "cat3Contraction"],
  cor_intercept_exp = cor_re["(Intercept)", "cat3Expansion"],
  cor_con_exp = cor_re["cat3Contraction", "cat3Expansion"],
  expected_pct_species_reverse = 100 * pnorm(g[["est"]] / sd_ce),
  n_species_con_and_exp = nrow(both),
  pct_species_fit_con_lt_exp = 100 * mean(both$fit_con_minus_exp < 0),
  resid_pct_beyond_1pct = pct_beyond_1pct, resid_pct_beyond_0.1pct = pct_beyond_0.1pct,
  r2_within_species = r2_within, r2_species = r2_species,
  group_lrt_chisq = lrt$Chisq[2], group_lrt_df = lrt$`Chi Df`[2], group_lrt_p = lrt$`Pr(>Chisq)`[2],
  group_delta_aic = lrt$AIC[1] - lrt$AIC[2],
  group_sd_con_shared = sd_shared[["cat3Contraction"]], group_sd_con_by_group = sd_group[["cat3Contraction"]],
  group_sd_exp_shared = sd_shared[["cat3Expansion"]], group_sd_exp_by_group = sd_group[["cat3Expansion"]]
)

write.csv(summary_out, out_summary, row.names = FALSE)
write.csv(species_out, out_species, row.names = FALSE)
write.csv(group_out, out_group, row.names = FALSE)
cat("\nWrote:\n ", out_summary, "\n ", out_species, "\n ", out_group, "\n")
cat("\nDone.\n")
