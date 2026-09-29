# 4f_sdm_category_mixed_models.R
#
# Tests whether SDM habitat-change category (Contraction / Stable /
# Expansion) predicts route-level trend WITHIN species, using linear mixed
# models (lme4::lmer) with species as the unit of replication, plus a
# robust Student-t version of the key models (glmmTMB).
#
# Why this exists alongside 4c/4d/4e: 4c's Kruskal-Wallis / Wilcoxon tests
# pool every route of every species as if independent, and 4d's category
# means mix two things -- the within-species effect the hypothesis is about,
# and species composition (species that are doing well everywhere may simply
# have more Expansion routes). Here each species gets its own baseline
# (random intercept) AND its own Contraction/Expansion effects (random
# slopes), so the fixed effects are the average WITHIN-species effect across
# species, and their SEs reflect how much that effect varies between species
# (~440 species), not 213k routes.
#
# The headline contrast throughout is CONTRACTION - EXPANSION (C - E, %/yr):
# the hypothesis predicts it is NEGATIVE (routes in a species' Contraction
# areas trend worse than its Expansion routes). "gap" columns below are all
# C - E. The two sensitivity-check slopes (weighted_score, grad7) are per
# step TOWARD Expansion, so they are predicted POSITIVE.
#
# Runs anthro x RCP4.5 by default (model_tags / scenarios under Settings);
# the output CSVs keep model/scenario columns so more combinations can be
# added there. Per model tag x scenario:
#   1. Three nested fits of trend ~ cat3 (reference level = Stable), to show
#      what each ingredient does to the C - E gap:
#        naive      lm(trend ~ cat3)                         routes independent
#        ri         lmer(trend ~ cat3 + (1 | species))       species baseline only
#        rs (MAIN)  lmer(trend ~ cat3 + (1 + cat3 | species)) baseline + category effects vary by species
#      C - E = beta_C - beta_E, SE = sqrt(L V L'), L = (0, 1, -1).
#      p-values use the normal approximation (t as z) -- lme4 gives none, and
#      with ~440 species behind the slope terms the approximation is safe.
#   2. Sensitivity checks on the main model:
#        weighted   ordinal score (-1/0/1) with inverse-variance weights from
#                   each route's 90% CI (se = (uci - lci) / 3.29, floored at
#                   0.25, weights normalised to mean 1)
#        grad7      the full 7-level raster code (1..7) as a numeric gradient
#        sign test  species with >=3 Contraction and >=3 Expansion routes:
#                   share with raw mean(C) < mean(E), binom.test vs 50%
#        Student-t  the main model refit in glmmTMB with t_family() errors
#                   (df estimated). Route trends are heavy-tailed (residual
#                   excess kurtosis ~6), which the Gaussian model ignores;
#                   AIC compares it against the ML refit of the main model.
#                   Also: its random-effect SDs/correlations, a tail check
#                   (share of residuals beyond the t(df) 1% / 0.1% cutoffs,
#                   vs the same check for the Gaussian fit against the
#                   normal), and the weighted / grad7 checks refit with
#                   t errors.
#   3. Predictive strength: R2 of cat3 on species-demeaned trend (within-
#      species R2), vs R2 of species identity alone.
#   4. Bird group, two ways (arctic dropped: 1 species, 2 Stable routes, no
#      category contrast possible):
#        Option A   group REPLACES species: (1 + cat3 | group). Kept only to
#                   show why it's a poor design (routes of a species treated
#                   as independent again; 11 groups too few for a variance).
#        Option B   species within groups (recommended):
#                   trend ~ 0 + group + group:cat3 + (1 + cat3 | species)
#                   per-group C - E with BH q across groups, and an ML
#                   likelihood-ratio test vs a shared-category-effect model
#                   (trend ~ 0 + group + cat3 + (1 + cat3 | species)) --
#                   ML, not REML, because REML likelihoods can't compare
#                   models with different fixed effects. Fit both Gaussian
#                   (lmer) and Student-t (glmmTMB); Option A too.
#
# Write-ups of the anthro / RCP4.5 results (hand-built HTML, not generated
# by this script):
#   output/reports/sdm_category_trend_report_anthro_rcp45_2010_2025.html
#     Gaussian LMM as the main model, Student-t as a robustness check.
#   output/reports/sdm_category_trend_report_student_t_anthro_rcp45_2010_2025.html
#     Student-t model as the main model, Gaussian LMM as a sensitivity check.
#
# Reads:  output/species_routes_covariates/per_species_sdm/*_route_trends_sdm.csv
#         (written by 3c_add_SDM_covariates.R)
# Writes: output/files/sdm_category_mixed_model_summary_2010_2025.csv
#           One row per model x scenario: data counts, the three C - E gaps,
#           main-model fixed effects + random-effect SDs, residual tail
#           diagnostics, Student-t gap/df/AIC, sensitivity checks, R2s, sign
#           test, Option A gap, Option B likelihood-ratio tests.
#         output/files/sdm_category_mixed_model_species_effects_2010_2025.csv
#           One row per model x scenario x species: route counts per category,
#           raw category means, and the main model's species-level fitted
#           Contraction / Stable / Expansion trends and C - E effect
#           (fixed + random, i.e. coef(fit); shrunk toward the average), plus
#           the Student-t model's species-level C - E.
#         output/files/sdm_category_mixed_model_by_group_2010_2025.csv
#           One row per model x scenario x group: Option B per-group
#           Stable baseline, C - S, E - S, C - E (SE, 95% CI, p, BH q) from
#           the Gaussian fit, the same C - E columns from the Student-t fit
#           (t_ prefix), and Option A's group-level C - E for comparison.

library(here)
library(tidyverse)
library(lme4)
library(glmmTMB)

here::i_am("4f_sdm_category_mixed_models.R")

# Settings -------------------------------------------------------------------
model_tags <- c("anthro")   # add "base" to also run the no-covariate trend model
scenarios  <- c("rcp45")    # add "rcp85" to also run the RCP8.5 categories
first_year <- 2010
last_year  <- 2025
require_route_converged <- TRUE   # same rule as 4c
min_routes_sign_test <- 3         # per category, for the species sign test
se_floor <- 0.25                  # floor on route trend SE for the weighted check
excluded_groups <- c("arctic")    # no Contraction/Expansion routes at all

in_dir  <- here::here("output", "species_routes_covariates", "per_species_sdm")
out_dir <- here::here("output", "files")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
suffix <- paste0("_", first_year, "_", last_year, ".csv")
out_summary <- file.path(out_dir, paste0("sdm_category_mixed_model_summary", suffix))
out_species <- file.path(out_dir, paste0("sdm_category_mixed_model_species_effects", suffix))
out_group   <- file.path(out_dir, paste0("sdm_category_mixed_model_by_group", suffix))

ctl <- lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
p_norm <- function(z) 2 * pnorm(-abs(z))
z975 <- qnorm(0.975)

# Linear contrast of fixed effects: estimate and SE.
contrast <- function(b, V, l) {
  l <- l[names(b)]
  l[is.na(l)] <- 0
  c(est = sum(l * b), se = sqrt(drop(t(l) %*% V %*% l)))
}
# Weights for Contraction - Expansion.
ce_weights <- function(con_name, exp_name) {
  setNames(c(1, -1), c(con_name, exp_name))
}
# glmmTMB's fixef()/vcov() return lists; take the conditional model.
tmb_b <- function(f) fixef(f)$cond
tmb_V <- function(f) as.matrix(vcov(f)$cond)
tmb_t_df <- function(f) {
  tryCatch(unname(family_params(f)[1]), error = function(e) exp(unname(f$fit$par["psi"])))
}

# Per-group C - E table from a group:cat3 fit (Gaussian or Student-t).
group_ce_table <- function(b, V, groups) {
  map_dfr(groups, function(g) {
    nm_s <- paste0("group", g)
    nm_c <- paste0("group", g, ":cat3Contraction")
    nm_e <- paste0("group", g, ":cat3Expansion")
    if (!all(c(nm_c, nm_e) %in% names(b))) return(NULL)  # group lacks a category
    ce <- contrast(b, V, ce_weights(nm_c, nm_e))
    tibble(group = g, stable = b[[nm_s]],
           con_minus_stable = b[[nm_c]], con_minus_stable_se = sqrt(V[nm_c, nm_c]),
           exp_minus_stable = b[[nm_e]], exp_minus_stable_se = sqrt(V[nm_e, nm_e]),
           con_minus_exp = ce[["est"]], con_minus_exp_se = ce[["se"]])
  }) %>%
    mutate(ci_lo = con_minus_exp - z975 * con_minus_exp_se,
           ci_hi = con_minus_exp + z975 * con_minus_exp_se,
           p = p_norm(con_minus_exp / con_minus_exp_se),
           q_bh = p.adjust(p, method = "BH"))
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

summary_rows <- list()
species_rows <- list()
group_rows   <- list()

for (m in model_tags) {
  for (sc in scenarios) {
    lab <- paste(m, sc)
    cat("\n==================", lab, "==================\n")

    x <- all_routes %>%
      filter(model == m) %>%
      transmute(species, species_code, group, trend, trend_lci, trend_uci, cat7 = .data[[sc]]) %>%
      filter(!is.na(cat7), cat7 %in% 1:7) %>%
      mutate(cat3 = factor(case_when(cat7 <= 3 ~ "Contraction", cat7 == 4 ~ "Stable", TRUE ~ "Expansion"),
                           levels = c("Stable", "Contraction", "Expansion")),
             score = case_when(cat3 == "Contraction" ~ -1, cat3 == "Stable" ~ 0, TRUE ~ 1),
             se = (trend_uci - trend_lci) / (2 * qnorm(0.95)))
    if (nrow(x) == 0) { cat("No rows -- skipped\n"); next }
    n_cat_by_sp <- x %>% group_by(species_code) %>% summarise(n_cat = n_distinct(cat3), .groups = "drop")
    cat("routes:", nrow(x), " species:", n_distinct(x$species_code),
        " species with >=2 categories:", sum(n_cat_by_sp$n_cat >= 2), "\n")

    # 1. Naive / random intercept / random intercept + slopes ---------------
    f_naive <- lm(trend ~ cat3, data = x)
    f_ri    <- lmer(trend ~ cat3 + (1 | species_code), data = x, control = ctl)
    f_rs    <- lmer(trend ~ cat3 + (1 + cat3 | species_code), data = x, control = ctl)
    if (isSingular(f_rs)) warning(lab, ": main model fit is singular")

    l_ce <- ce_weights("cat3Contraction", "cat3Expansion")
    g_naive <- contrast(coef(f_naive), vcov(f_naive), l_ce)
    g_ri    <- contrast(fixef(f_ri), as.matrix(vcov(f_ri)), l_ce)
    g_rs    <- contrast(fixef(f_rs), as.matrix(vcov(f_rs)), l_ce)
    fe <- coef(summary(f_rs))
    vc <- as.data.frame(VarCorr(f_rs))
    sd_of <- function(v1) vc$sdcor[vc$grp == "species_code" & vc$var1 == v1 & is.na(vc$var2)]
    cov_ce <- vc$vcov[vc$grp == "species_code" & vc$var1 == "cat3Contraction" & vc$var2 %in% "cat3Expansion"]
    sd_ce <- sqrt(sd_of("cat3Contraction")^2 + sd_of("cat3Expansion")^2 - 2 * cov_ce)
    cat(sprintf("C - E gap: naive %.2f | random intercept %.2f | main model %.2f (SE %.2f, p = %.2g)\n",
                g_naive["est"], g_ri["est"], g_rs["est"], g_rs["se"], p_norm(g_rs["est"] / g_rs["se"])))

    # Residual tails of the Gaussian main model: share of standardised
    # residuals beyond the cutoffs a normal distribution exceeds 1% / 0.1% of
    # the time (two-sided).
    r_std <- residuals(f_rs) / sigma(f_rs)
    resid_excess_kurtosis <- mean(r_std^4) / mean(r_std^2)^2 - 3
    resid_pct_beyond_1pct   <- 100 * mean(abs(r_std) > qnorm(0.995))
    resid_pct_beyond_0.1pct <- 100 * mean(abs(r_std) > qnorm(0.9995))

    # 2. Sensitivity checks -------------------------------------------------
    x <- x %>% mutate(w = 1 / pmax(se, se_floor)^2, w = w / mean(w))
    f_wtd   <- lmer(trend ~ score + (1 + score | species_code), data = x, weights = w, control = ctl)
    f_grad7 <- lmer(trend ~ cat7 + (1 + cat7 | species_code), data = x, control = ctl)

    # Student-t errors (robust to the heavy-tailed route trends)
    f_t <- glmmTMB(trend ~ cat3 + (1 + cat3 | species_code), data = x, family = t_family())
    g_t <- contrast(tmb_b(f_t), tmb_V(f_t), l_ce)
    fe_t <- coef(summary(f_t))$cond
    aic_gauss_ml <- AIC(refitML(f_rs))
    t_df_hat <- tmb_t_df(f_t)
    cat(sprintf("Student-t: C - E %.2f (SE %.2f, p = %.2g), df = %.2f, AIC lower than Gaussian by %.0f\n",
                g_t["est"], g_t["se"], p_norm(g_t["est"] / g_t["se"]), t_df_hat, aic_gauss_ml - AIC(f_t)))
    # Student-t random effects: SDs, correlations, derived C - E spread
    vc_t <- VarCorr(f_t)$cond$species_code
    sd_t <- attr(vc_t, "stddev")
    cor_t <- attr(vc_t, "correlation")
    sd_ce_t <- sqrt(vc_t["cat3Contraction", "cat3Contraction"] + vc_t["cat3Expansion", "cat3Expansion"] -
                      2 * vc_t["cat3Contraction", "cat3Expansion"])
    # Student-t residual tails: share beyond the cutoffs a t(df_hat)
    # distribution exceeds 1% / 0.1% of the time -- checks that the fitted t
    # actually captures the tails, not just that it beats the normal.
    r_t <- residuals(f_t) / sigma(f_t)
    t_resid_pct_beyond_1pct   <- 100 * mean(abs(r_t) > qt(0.995, t_df_hat))
    t_resid_pct_beyond_0.1pct <- 100 * mean(abs(r_t) > qt(0.9995, t_df_hat))
    # Sensitivity checks refit with Student-t errors
    f_wtd_t   <- glmmTMB(trend ~ score + (1 + score | species_code), data = x, weights = w, family = t_family())
    f_grad7_t <- glmmTMB(trend ~ cat7 + (1 + cat7 | species_code), data = x, family = t_family())

    sp_raw <- x %>%
      group_by(species_code) %>%
      summarise(species = first(species), group = first(group), n_routes = n(),
                n_con = sum(cat3 == "Contraction"), n_stable = sum(cat3 == "Stable"),
                n_exp = sum(cat3 == "Expansion"),
                raw_con = if (n_con > 0) mean(trend[cat3 == "Contraction"]) else NA_real_,
                raw_stable = if (n_stable > 0) mean(trend[cat3 == "Stable"]) else NA_real_,
                raw_exp = if (n_exp > 0) mean(trend[cat3 == "Expansion"]) else NA_real_,
                .groups = "drop")
    sign_sp <- sp_raw %>% filter(n_con >= min_routes_sign_test, n_exp >= min_routes_sign_test) %>%
      mutate(diff = raw_con - raw_exp)
    bt <- binom.test(sum(sign_sp$diff < 0), nrow(sign_sp))

    # 3. Predictive strength ------------------------------------------------
    x_dm <- x %>% left_join(n_cat_by_sp, by = "species_code") %>% filter(n_cat >= 2) %>%
      group_by(species_code) %>% mutate(trend_dm = trend - mean(trend)) %>% ungroup()
    r2_within  <- summary(lm(trend_dm ~ cat3, data = x_dm))$r.squared
    r2_species <- summary(lm(trend ~ species_code, data = x))$r.squared

    # Species-level estimates (Gaussian main model + Student-t) -------------
    sp_coef <- coef(f_rs)$species_code %>%
      rownames_to_column("species_code") %>%
      as_tibble() %>%
      transmute(species_code,
                fit_stable = `(Intercept)`,
                fit_con = `(Intercept)` + cat3Contraction,
                fit_exp = `(Intercept)` + cat3Expansion,
                fit_con_minus_exp = cat3Contraction - cat3Expansion)
    sp_coef_t <- coef(f_t)$cond$species_code %>%
      rownames_to_column("species_code") %>%
      as_tibble() %>%
      transmute(species_code,
                t_fit_stable = `(Intercept)`,
                t_fit_con = `(Intercept)` + cat3Contraction,
                t_fit_exp = `(Intercept)` + cat3Expansion,
                t_fit_con_minus_exp = cat3Contraction - cat3Expansion)
    sp_out <- sp_raw %>%
      left_join(sp_coef, by = "species_code") %>%
      left_join(sp_coef_t, by = "species_code") %>%
      mutate(model = m, scenario = sc, .before = 1)
    species_rows[[lab]] <- sp_out
    both <- sp_out %>% filter(n_con >= 1, n_exp >= 1)

    # 4. Bird group ---------------------------------------------------------
    xg <- x %>% filter(!group %in% excluded_groups) %>% mutate(group = factor(group))

    # Option A: group replaces species
    f_a <- lmer(trend ~ cat3 + (1 + cat3 | group), data = xg, control = ctl)
    g_a <- contrast(fixef(f_a), as.matrix(vcov(f_a)), l_ce)
    a_groups <- coef(f_a)$group %>% rownames_to_column("group") %>% as_tibble() %>%
      transmute(group, optA_con_minus_exp = cat3Contraction - cat3Expansion)
    f_a_t <- glmmTMB(trend ~ cat3 + (1 + cat3 | group), data = xg, family = t_family())
    g_a_t <- contrast(tmb_b(f_a_t), tmb_V(f_a_t), l_ce)
    fe_a_t <- coef(summary(f_a_t))$cond

    # Option B: species within groups, group-specific category effects
    f_b  <- lmer(trend ~ 0 + group + group:cat3 + (1 + cat3 | species_code), data = xg, control = ctl, REML = FALSE)
    f_b0 <- lmer(trend ~ 0 + group + cat3 + (1 + cat3 | species_code), data = xg, control = ctl, REML = FALSE)
    lrt <- anova(f_b0, f_b)
    grp_tab <- group_ce_table(fixef(f_b), as.matrix(vcov(f_b)), levels(xg$group))

    # Option B with Student-t errors (glmmTMB fits by ML, so the LRT is direct)
    f_bt  <- glmmTMB(trend ~ 0 + group + group:cat3 + (1 + cat3 | species_code), data = xg, family = t_family())
    f_bt0 <- glmmTMB(trend ~ 0 + group + cat3 + (1 + cat3 | species_code), data = xg, family = t_family())
    lrt_t <- anova(f_bt0, f_bt)
    grp_tab_t <- group_ce_table(tmb_b(f_bt), tmb_V(f_bt), levels(xg$group)) %>%
      select(group, con_minus_exp, con_minus_exp_se, ci_lo, ci_hi, p, q_bh) %>%
      rename_with(~ paste0("t_", .x), -group)
    # How much of the between-species spread in category effects do groups
    # explain? Species-level SDs with a shared vs group-specific category effect.
    sd_bt0 <- attr(VarCorr(f_bt0)$cond$species_code, "stddev")
    sd_bt  <- attr(VarCorr(f_bt)$cond$species_code, "stddev")

    grp_counts <- xg %>% group_by(group) %>%
      summarise(n_species = n_distinct(species_code), n_routes = n(),
                n_con = sum(cat3 == "Contraction"), n_stable = sum(cat3 == "Stable"),
                n_exp = sum(cat3 == "Expansion"), .groups = "drop") %>%
      mutate(group = as.character(group))
    group_rows[[lab]] <- grp_tab %>%
      left_join(grp_tab_t, by = "group") %>%
      left_join(grp_counts, by = "group") %>%
      left_join(a_groups, by = "group") %>%
      mutate(model = m, scenario = sc, .before = 1) %>%
      arrange(con_minus_exp)
    cat(sprintf("Group differences (Option B LRT): Gaussian chisq = %.1f, df = %d, p = %.2g | Student-t chisq = %.1f, p = %.2g\n",
                lrt$Chisq[2], lrt$Df[2], lrt$`Pr(>Chisq)`[2], lrt_t$Chisq[2], lrt_t$`Pr(>Chisq)`[2]))
    print(as.data.frame(group_rows[[lab]] %>%
                          select(group, n_species, con_minus_exp, ci_lo, ci_hi, q_bh,
                                 t_con_minus_exp, t_ci_lo, t_ci_hi, t_q_bh)),
          digits = 3, row.names = FALSE)

    # Summary row -----------------------------------------------------------
    summary_rows[[lab]] <- tibble(
      model = m, scenario = sc,
      n_routes = nrow(x), n_species = n_distinct(x$species_code),
      n_species_2cat = sum(n_cat_by_sp$n_cat >= 2),
      n_con = sum(x$cat3 == "Contraction"), n_stable = sum(x$cat3 == "Stable"), n_exp = sum(x$cat3 == "Expansion"),
      gap_naive = g_naive[["est"]], gap_naive_se = g_naive[["se"]],
      gap_ri = g_ri[["est"]], gap_ri_se = g_ri[["se"]],
      gap_main = g_rs[["est"]], gap_main_se = g_rs[["se"]],
      gap_main_ci_lo = g_rs[["est"]] - z975 * g_rs[["se"]],
      gap_main_ci_hi = g_rs[["est"]] + z975 * g_rs[["se"]],
      gap_main_p = p_norm(g_rs[["est"]] / g_rs[["se"]]),
      main_singular = isSingular(f_rs),
      stable_avg_species = fe["(Intercept)", 1], stable_avg_species_se = fe["(Intercept)", 2],
      con_minus_stable = fe["cat3Contraction", 1], con_minus_stable_se = fe["cat3Contraction", 2],
      exp_minus_stable = fe["cat3Expansion", 1], exp_minus_stable_se = fe["cat3Expansion", 2],
      sd_species_intercept = sd_of("(Intercept)"), sd_species_con = sd_of("cat3Contraction"),
      sd_species_exp = sd_of("cat3Expansion"), sd_species_con_minus_exp = sd_ce,
      sd_residual = sigma(f_rs),
      resid_excess_kurtosis = resid_excess_kurtosis,
      resid_pct_beyond_3sd = 100 * mean(abs(r_std) > 3),
      resid_pct_beyond_1pct = resid_pct_beyond_1pct, resid_pct_beyond_0.1pct = resid_pct_beyond_0.1pct,
      expected_pct_species_reverse = 100 * pnorm(g_rs[["est"]] / sd_ce),
      n_species_con_and_exp = nrow(both),
      pct_species_fit_con_lt_exp = 100 * mean(both$fit_con_minus_exp < 0),
      t_gap = g_t[["est"]], t_gap_se = g_t[["se"]],
      t_gap_ci_lo = g_t[["est"]] - z975 * g_t[["se"]],
      t_gap_ci_hi = g_t[["est"]] + z975 * g_t[["se"]],
      t_gap_p = p_norm(g_t[["est"]] / g_t[["se"]]),
      t_con_minus_stable = fe_t["cat3Contraction", 1], t_con_minus_stable_se = fe_t["cat3Contraction", 2],
      t_exp_minus_stable = fe_t["cat3Expansion", 1], t_exp_minus_stable_se = fe_t["cat3Expansion", 2],
      t_stable_avg_species = fe_t["(Intercept)", 1], t_stable_avg_species_se = fe_t["(Intercept)", 2],
      t_df = t_df_hat, t_sigma = sigma(f_t),
      t_converged = f_t$fit$convergence == 0 && isTRUE(f_t$sdr$pdHess),
      t_sd_species_intercept = sd_t[["(Intercept)"]], t_sd_species_con = sd_t[["cat3Contraction"]],
      t_sd_species_exp = sd_t[["cat3Expansion"]], t_sd_species_con_minus_exp = sd_ce_t,
      t_cor_intercept_con = cor_t["(Intercept)", "cat3Contraction"],
      t_cor_intercept_exp = cor_t["(Intercept)", "cat3Expansion"],
      t_cor_con_exp = cor_t["cat3Contraction", "cat3Expansion"],
      t_expected_pct_species_reverse = 100 * pnorm(g_t[["est"]] / sd_ce_t),
      t_resid_pct_beyond_1pct = t_resid_pct_beyond_1pct, t_resid_pct_beyond_0.1pct = t_resid_pct_beyond_0.1pct,
      t_pct_species_fit_con_lt_exp = 100 * mean(both$t_fit_con_minus_exp < 0),
      t_species_effect_cor_with_gaussian = cor(both$fit_con_minus_exp, both$t_fit_con_minus_exp),
      aic_gaussian_ml = aic_gauss_ml, aic_t = AIC(f_t),
      t_weighted_score = tmb_b(f_wtd_t)[["score"]], t_weighted_score_se = sqrt(tmb_V(f_wtd_t)["score", "score"]),
      t_grad7 = tmb_b(f_grad7_t)[["cat7"]], t_grad7_se = sqrt(tmb_V(f_grad7_t)["cat7", "cat7"]),
      weighted_score = fixef(f_wtd)[["score"]], weighted_score_se = coef(summary(f_wtd))["score", 2],
      grad7 = fixef(f_grad7)[["cat7"]], grad7_se = coef(summary(f_grad7))["cat7", 2],
      r2_naive = summary(f_naive)$r.squared, r2_species = r2_species, r2_within_species = r2_within,
      sign_test_n_species = nrow(sign_sp), sign_test_pct_con_lt_exp = 100 * mean(sign_sp$diff < 0),
      sign_test_p = bt$p.value, sign_test_median_diff = median(sign_sp$diff),
      optA_gap = g_a[["est"]], optA_gap_se = g_a[["se"]], optA_singular = isSingular(f_a),
      optB_lrt_chisq = lrt$Chisq[2], optB_lrt_df = lrt$Df[2], optB_lrt_p = lrt$`Pr(>Chisq)`[2],
      optB_delta_aic = lrt$AIC[1] - lrt$AIC[2],
      # glmmTMB's anova() puts the parameter count in "Df" and the test df in "Chi Df"
      optB_t_lrt_chisq = lrt_t$Chisq[2], optB_t_lrt_df = lrt_t$`Chi Df`[2], optB_t_lrt_p = lrt_t$`Pr(>Chisq)`[2],
      optB_t_delta_aic = lrt_t$AIC[1] - lrt_t$AIC[2],
      optA_t_gap = g_a_t[["est"]], optA_t_gap_se = g_a_t[["se"]],
      optA_t_con_minus_stable = fe_a_t["cat3Contraction", 1], optA_t_con_minus_stable_se = fe_a_t["cat3Contraction", 2],
      optA_t_exp_minus_stable = fe_a_t["cat3Expansion", 1], optA_t_exp_minus_stable_se = fe_a_t["cat3Expansion", 2],
      optA_t_sigma = sigma(f_a_t),
      optB_t_sd_con_shared = sd_bt0[["cat3Contraction"]], optB_t_sd_con_by_group = sd_bt[["cat3Contraction"]],
      optB_t_sd_exp_shared = sd_bt0[["cat3Expansion"]], optB_t_sd_exp_by_group = sd_bt[["cat3Expansion"]]
    )
  }
}

summary_out <- bind_rows(summary_rows)
write.csv(summary_out, out_summary, row.names = FALSE)
write.csv(bind_rows(species_rows), out_species, row.names = FALSE)
write.csv(bind_rows(group_rows), out_group, row.names = FALSE)
cat("\nWrote:\n ", out_summary, "\n ", out_species, "\n ", out_group, "\n")

cat("\n=== Contraction - Expansion gap (%/yr), per model x scenario ===\n")
print(as.data.frame(summary_out %>%
  select(model, scenario, gap_naive, gap_ri, gap_main, gap_main_ci_lo, gap_main_ci_hi,
         t_gap, t_gap_ci_lo, t_gap_ci_hi, t_df, r2_within_species, sign_test_pct_con_lt_exp,
         optB_lrt_p, optB_t_lrt_p)),
  digits = 3, row.names = FALSE)

cat("\nDone.\n")
