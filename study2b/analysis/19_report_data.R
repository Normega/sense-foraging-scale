# 19: the data behind senseforaging.com/research/study2b/results. Builds the site's data:
# aggregates only, read from the real-mode outputs of 01-18, 20 and 21 (and the first-700
# sensitivity run's 09_decisions.csv). Writes study2bResults.json to SF2B_REPORT_OUT (default:
# the Handbook checkout's src/pages/research/). Fits nothing; rerun after any of them changes.
#
# Status of each block, as the page labels it:
#   registered   01-09 (H1-H6, the registered exploratory H8 and descriptives in 08)
#   deviation    10 (decision 46's Brief MAIA-2 repair) and the first-700 run (decision 45)
#   exploratory  11-18, 20 and 21: specified after the confirmatory run, not in osf.io/d4a5u
#
# A few facts are not in any output file: they are the decision log's record of data
# collection (study2/decisions_after_2a.md, decisions 44-46) and are written under
# `recorded`, each with its source.
#
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/19_report_data.R

# Set Up ---------
## Load libraries ---------
packages <- c("jsonlite", "data.table")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

## Settings ---------
source("study2b/analysis/config.R")
stopifnot(SF2B$mode == "real", SF2B$first_n == 0)
dest <- Sys.getenv("SF2B_REPORT_OUT", "F:/gits/Handbook/src/pages/research/study2bResults.json")

o  <- function(f) fread(sf2b_out(f))
oj <- function(f) fromJSON(sf2b_out(f), simplifyVector = FALSE)
# Values keep 4 decimals (or the source's fewer), so the page's own rounding is not a double rounding.
r  <- function(x, d = 4) if (is.null(x) || length(x) == 0) NULL else round(as.numeric(x), d)
rows <- function(dt) lapply(seq_len(nrow(dt)), function(i) as.list(dt[i]))
# The one row of dt where cond holds (cond is evaluated in dt, then in the caller's frame).
one <- function(dt, cond) {
  keep <- eval(substitute(cond), dt, parent.frame())
  x <- dt[which(keep)]
  if (nrow(x) != 1) stop("expected one row, got ", nrow(x), ": ", deparse(substitute(cond)))
  x
}

# Short measure names for the page (names only, never item wording).
LAB <- c(
  ffmq_observing = "FFMQ Observing", ffmq_describing = "FFMQ Describing", ffmq_acting = "FFMQ Acting with Awareness",
  ffmq_nonjudging = "FFMQ Non-judging", ffmq_nonreactivity = "FFMQ Non-reactivity",
  sbi_moment = "SBI Savoring the Moment", bfi_aesthetic = "BFI-2 Aesthetic Sensitivity",
  bfi_intellectual = "BFI-2 Intellectual Curiosity", bfi_creative = "BFI-2 Creative Imagination",
  bfi_openmindedness = "BFI-2 Open-Mindedness (total)", eq_decentering = "EQ Decentering",
  samd_md = "SAMD Meaning-Detecting", bsss8 = "BSSS-8 Sensation Seeking", ncs6 = "NCS-6 Need for Cognition",
  spsq_sd = "SPSQ Sensory Discomfort", asi_pc = "ASI-3 Physical Concerns",
  bidr_sde = "BIDR Self-Deceptive Enhancement", bidr_im = "BIDR Impression Management",
  phq4_total = "PHQ-4", pss4 = "PSS-4 stress", lifesat = "Life satisfaction",
  crit_pause_days = "Past-week pause days", engaged_pause = "Engaged pause (in session)",
  crit_notice_first = "Past week: noticed first", crit_outdoors = "Past week: time outdoors",
  maia_noticing = "Brief MAIA-2 Noticing", maia_notdistracting = "Brief MAIA-2 Not-Distracting",
  maia_notworrying = "Brief MAIA-2 Not-Worrying", maia_attnregulation = "Brief MAIA-2 Attention Regulation",
  maia_emoawareness = "Brief MAIA-2 Emotional Awareness", maia_selfregulation = "Brief MAIA-2 Self-Regulation",
  maia_bodylistening = "Brief MAIA-2 Body Listening", maia_trusting = "Brief MAIA-2 Trusting",
  vig_self = "Self-rating of usual mode (vignettes)", separation = "Vignette separation",
  sf8 = "The 8", grp = "Practice group (current vs never)")
lab <- function(x) unname(ifelse(x %in% names(LAB), LAB[x], x))

# Sample and flow (01, 02) ---------
flow <- o("01_flow.csv"); flow[, step := trimws(step)]
fl <- function(s) flow[step == s, n]
flags <- o("01_flags.csv")[analysable == TRUE]
pc <- o("01_pause_check.csv")
smp <- oj("02_sample.json")
sample <- list(
  enrolled = fl("enrollments"), staff = fl("excluded: staff test enrollment"),
  no_scale = fl("excluded: no submitted scale"), failed_checks = fl("excluded: failed two or more attention checks"),
  analysable = fl("analysable (H1)"), reached_validity = fl("reached the first half of the validity block"),
  reached_pause = fl("reached the pause offer (H4b)"), reached_end = fl("reached the end of the session"),
  flagged = flags$any_flag, one_check = flags$one_check, speeder = flags$speeder, longstring = flags$longstring,
  pna = flags$pna_over_2,
  age_mean = smp$age_mean, age_sd = smp$age_sd,
  gender = smp$gender[c("man", "woman", "non_binary")],
  practice = smp$practice_group[c("current", "past_or_occasional", "never")], dedicated = smp$dedicated,
  tells_apart = smp$tells_apart,
  pause = list(offered = pc$reached_offer, took = pc$took_pause, engaged = pc$engaged, note = pc$wrote_a_note,
               cap = pc$ended_by_cap, nudged = pc$nudged),
  device = smp$pause_device[c("desktop", "phone", "tablet")]
)

recorded <- list(
  source = "study2/decisions_after_2a.md, decisions 44-46 (the decision log; not an analysis output)",
  pilot_date = "2026-11-13", main_date = "2026-11-15", filed = "2026-11-14", prereg = "osf.io/d4a5u", prereg_commit = "d3e60c9",
  completed = 723, pilot_completed = 61, main_completed = 662, reentry = 5,
  pilot_engaged = 2, h4b_power_at_pilot_pct = 49, prolific_max_min = 48, prolific_max_registered_min = 90
)

# Confirmatory decisions (09) and the first-700 sensitivity run ---------
dec <- o("09_decisions.csv")
d700 <- fread(sf2b_path("output", "real_first700", "09_decisions.csv"))
stopifnot(identical(dec$test, d700$test))
first_num <- function(s) suppressWarnings(as.numeric(regmatches(s, regexpr("-?[0-9]+\\.[0-9]+", s))))
num_full <- vapply(dec$estimate, function(s) { x <- first_num(s); if (length(x)) x else NA_real_ }, 0)
num_700  <- vapply(d700$estimate, function(s) { x <- first_num(s); if (length(x)) x else NA_real_ }, 0)
f700 <- fromJSON(sf2b_path("output", "real_first700", "03_h1_decision.json"), simplifyVector = FALSE)
s700 <- fread(sf2b_path("output", "real_first700", "01_flow.csv"))
decisions <- list(
  rows = lapply(seq_len(nrow(dec)), function(i) list(hypothesis = dec$hypothesis[i], test = dec$test[i], result = dec$result[i],
                                                       estimate = dec$estimate[i], first700 = d700$result[i])),
  first700_same = all(dec$result == d700$result),
  first700_max_diff = r(max(abs(num_full - num_700), na.rm = TRUE)),
  first700_analysable = f700$h1$n,
  first700_enrolled = s700[trimws(step) == "enrollments", n]
)

# H1 (03) ---------
h1j <- oj("03_h1_decision.json"); fit <- o("03_fit.csv"); ld <- o("03_loadings.csv")
fitrow <- function(m) { x <- one(fit, model == m); list(cfi = r(x$cfi.scaled), rmsea = r(x$rmsea.scaled), lo = r(x$rmsea.ci.lower.scaled),
                                                         hi = r(x$rmsea.ci.upper.scaled), srmr = r(x$srmr)) }
h <- h1j$h1; t10 <- h1j$the_10
lat <- o("04_latent_validity.csv")[sample == "all"]
l8 <- lat[form == "sf8"]; l10 <- lat[form == "sf10"]
m10 <- merge(l8[, .(score, r8 = est)], l10[, .(score, r10 = est)], by = "score")
h1 <- list(
  fit = c(fitrow("the 8, with agreement (H1)"), list(chisq = r(h$chisq, 1), df = h$df, max_residual = h$max_residual,
          freed_omega = h$freeing_it_changes_omega_by, freed_loading = h$freeing_it_changes_a_loading_by_up_to,
          usrmr = h$usrmr, closefit_p = h$closefit_p, omega = h$omega, alpha = h$alpha,
          agreement_share = h$agreement_share_of_item_variance, supported = h$supported)),
  without_agreement = fitrow("the 8, without agreement"),
  unflagged = c(fitrow("the 8, without flagged participants"), list(n = h1j$h1_unflagged$n, supported = h1j$h1_unflagged$supported)),
  the10 = c(fitrow("the 10, with agreement"), list(omega = t10$omega, alpha = t10$alpha, supported = t10$supported,
            max_validity_diff = r(max(abs(m10$r10 - m10$r8))),
            max_validity_diff_score = lab(m10$score[which.max(abs(m10$r10 - m10$r8))]))),
  loadings = lapply(seq_len(nrow(ld[form == 8])), function(i) { x <- ld[form == 8][i]
    list(item = x$item, code = x$source, reversed = x$content < 0, content = r(x$keyed_content), agreement = r(x$agreement),
         residual = r(x$residual_variance)) }),
  added_in_10 = lapply(c("A9", "sfc_24"), function(s) { x <- one(ld, form == 10 & source == s); list(code = s, content = r(x$keyed_content)) }),
  unscored = lapply(h1j$exploratory_statements, function(x) list(code = x$source, content = r(x$keyed_content), agreement = r(x$agreement))),
  agree_pct = rows(o("02_agreement_rates.csv")[, .(code = source, role, pct_agree)])
)
h1$loadings <- h1$loadings[order(-vapply(h1$loadings, `[[`, 0, "content"))]

# H2-H3 (04) and the MAIA repair (10) ---------
vd <- oj("04_validity_decisions.json")$sf8
mj <- oj("10_deviation_maia.json")$sf8
chart_groups <- list(
  list(g = "H2 neighbours", s = c("ffmq_observing", "maia_noticing", "sbi_moment", "bfi_aesthetic")),
  list(g = "Contrast partners and related", s = c("samd_md", "bfi_creative", "eq_decentering", "bfi_openmindedness", "ffmq_describing",
                                                  "bfi_intellectual", "ffmq_nonreactivity", "ffmq_acting", "ffmq_nonjudging")),
  list(g = "Discriminant (H3c, H3e)", s = c("bidr_sde", "ncs6", "bidr_im", "bsss8", "asi_pc", "spsq_sd")),
  list(g = "Outcomes and criteria", s = c("crit_pause_days", "lifesat", "engaged_pause", "phq4_total", "pss4")))
cell <- function(s) {
  if (s == "maia_noticing") { x <- mj$h2_noticing; return(list(score = s, label = lab(s), r = x$r, lo = x$lower, hi = x$upper, repair = TRUE)) }
  x <- one(l8, score == s); list(score = s, label = lab(s), r = r(x$est), lo = r(x$ci.lower), hi = r(x$ci.upper), repair = FALSE, prediction = x$prediction)
}
validity <- list(
  chart = lapply(chart_groups, function(g) list(group = g$g, rows = lapply(g$s, cell))),
  all = lapply(seq_len(nrow(l8)), function(i) list(score = l8$score[i], label = lab(l8$score[i]), r = r(l8$est[i]),
                                                     lo = r(l8$ci.lower[i]), hi = r(l8$ci.upper[i]), prediction = l8$prediction[i])),
  h2 = lapply(vd$h2, function(x) list(score = x$score, label = lab(x$score), r = x$r, lo = x$lower, hi = x$upper, h2 = x$h2, h3a = x$h3a, class = x$rc_class)),
  h2_holding = vd$H2_correlations_holding, H2 = vd$H2, H3a = vd$H3a,
  h3b = lapply(vd$h3b, function(x) list(contrast = x$contrast, first = lab(x$first), second = lab(x$second), r_first = x$r_first, r_second = x$r_second,
                                        diff = x$difference, lo = x$lower, hi = x$upper, holds = x$holds)),
  h3b_holding = vd$H3b_contrasts_holding, H3b = vd$H3b,
  h3c = lapply(vd$h3c, function(x) list(score = x$score, label = lab(x$score), r = x$r, lo = x$lower, hi = x$upper, holds = x$holds)), H3c = vd$H3c,
  h3d = list(R = vd$h3d$R, lo = vd$h3d$lower, hi = vd$h3d$upper, resamples = vd$h3d$resamples, resamples_ok = vd$h3d$resamples_ok,
             class = vd$h3d$rc_class, holds = vd$h3d$holds),
  h3e = lapply(vd$h3e, function(x) list(score = x$score, label = lab(x$score), r = x$r, lo = x$lower, hi = x$upper,
                                        bound_lo = x$bound_lower %||% NA, bound_hi = x$bound_upper, holds = x$holds)), H3e = vd$H3e,
  repair = list(noticing = mj$h2_noticing[c("r", "lower", "upper", "rc_class")], h2_holding = mj$H2$correlations_holding,
                H2 = mj$H2$holds, H3a = mj$H3a$holds, contrast4 = mj$contrast_4[c("r_first", "r_second", "difference", "lower", "upper", "holds")],
                h3b_holding = mj$H3b$contrasts_holding, H3b = mj$H3b$holds,
                std_loading_bound = r(abs(mj$maia_12$maia12_std_loading)), bound_active = mj$maia_12$bound_active)
)

# H4-H6 (05-07) ---------
h4j <- oj("05_h4_decision.json")$sf8
h5a <- o("06_h5a.csv")[form == "sf8"]; h5b <- one(o("06_h5b.csv"), form == "sf8")
kg <- o("07_known_groups.csv")[form == "sf8"]; ldf <- one(o("07_latent_difference.csv"), form == "sf8")
h4 <- list(a = lapply(h4j$a, r), b = lapply(h4j$b, r), b_logistic = list(or = r(h4j$b_logistic$odds_ratio_per_sd), lo = r(h4j$b_logistic$lower),
                                                                          hi = r(h4j$b_logistic$upper), covariate = h4j$b_logistic$covariate))
h4$a$p_holm <- signif(h4j$a$p_holm, 3); h4$a$supported <- h4j$a$supported; h4$b$supported <- h4j$b$supported
h5 <- list(a = rows(h5a[, .(outcome = lab(outcome), n, delta_r2 = r(delta_r2), p = r(p), p_holm = r(p_holm), supported)]),
           b = list(n = h5b$n, or = r(h5b$odds_ratio_per_sd), p = r(h5b$p), supported = h5b$supported))
h6 <- list(rows = rows(kg[, .(analysis, n_group, n_reference, mean_group = r(mean_group), mean_reference = r(mean_reference),
                              d = r(d), lo = r(d_lower), hi = r(d_upper), supported)]),
           latent = list(n = ldf$n, d = r(ldf$latent_d), lo = r(ldf$lower), hi = r(ldf$upper)))

# Vignettes, H8 and the registered exploratory analyses (08) ---------
ex <- oj("08_exploratory.json"); inv <- o("08_invariance.csv"); obsv <- o("04_observed_validity.csv")
vg <- ex$vignettes
h8 <- list(
  tells_apart = vg$tells_apart, share = r(vg$share), share_lo = r(vg$share_ci_wilson[[1]]), share_hi = r(vg$share_ci_wilson[[2]]),
  separation_mean = vg$separation_mean, separation_sd = r(vg$separation_sd),
  self_rating_r = r(vg$self_rating_with_sf8$r),
  means = lapply(names(vg$vignette_means), function(k) list(v = sub("vig_", "", k), kind = vg$vignette_kind[[k]], mean = vg$vignette_means[[k]])),
  invariance = rows(inv[, .(analysis, step, cfi = r(cfi), rmsea = r(rmsea), d_cfi = r(d_cfi, 4), noninvariant)]),
  groups = lapply(ex$invariance, function(x) list(groups = x$groups, last_step = x$last_step, released = length(unlist(x$released)))),
  country = list(analysed = ex$country_rule$analysed, six = ex$country_rule$groups$six_countries, other = ex$country_rule$groups$other, needed = 200),
  trial = list(r8 = r(ex$trial$sf8_with_trial_av10$r),
               observed = rows(obsv[score_of %in% c("trial_av10", "sf8") & with %in% c("ffmq_observing", "ncs6", "bidr_sde"),
                                    .(score_of, with = lab(with), r = r(r))])),
  bidr_dichotomous = list(sde = r(ex$bidr_dichotomous$sde$r), im = r(ex$bidr_dichotomous$im$r))
)

# The pause (11) ---------
pj <- oj("11_pause.json")
biv <- o("11_pause_bivariate.csv")
bv <- function(out, pred) one(biv, outcome == out & predictor == pred)
pr_t <- bv("took", "practice_current"); pr_e <- bv("engaged", "practice_current")
fq_t <- bv("took", "freq_code"); fq_e <- bv("engaged", "freq_code")
coef_of <- function(out, term) { x <- Filter(function(z) z$outcome == out && z$term == term, pj$model$coefficients)[[1]]; list(or = r(x$odds_ratio), lo = r(x$lower), hi = r(x$upper), p = signif(x$p, 3)) }
auc_of <- function(out) { x <- Filter(function(z) z$outcome == out, pj$model$auc)[[1]]; list(cv = r(x$auc_cv_mean), apparent = r(x$auc_apparent)) }
dur <- Filter(function(z) z$variable == "on_page_s", pj$paused$duration)[[1]]
pcor <- function(fam, w) { x <- Filter(function(z) z$family == fam && z$with == w, pj$paused$correlations)[[1]]; list(rho = x$rho, p_holm = x$p_holm) }
dose <- vapply(pj$measurement$dose, function(z) if (is.numeric(z$rho)) z$rho else NA_real_, 0)
inter <- o("11_pause_interactions.csv")
wave <- function(out) { x <- Filter(function(z) z$outcome == out, pj$practice_within_main_wave)[[1]]; list(n = x$n, or = r(x$odds_ratio), lo = r(x$lower), hi = r(x$upper)) }
pause <- list(
  offered = pj$base_rates$offered, took = pj$base_rates$took, engaged = pj$base_rates$engaged,
  took_rate = r(pj$base_rates$took_rate), took_lo = r(pj$base_rates$took_rate_ci_wilson[[1]]), took_hi = r(pj$base_rates$took_rate_ci_wilson[[2]]),
  engaged_rate = r(pj$base_rates$engaged_rate), engaged_lo = r(pj$base_rates$engaged_rate_ci_wilson[[1]]), engaged_hi = r(pj$base_rates$engaged_rate_ci_wilson[[2]]),
  tests_per_outcome = pj$bivariate$tests_per_outcome,
  practice = list(took = list(rate_current = r(pr_t$rate_if_1), rate_not = r(pr_t$rate_if_0), or = r(pr_t$odds_ratio), lo = r(pr_t$or_lower),
                              hi = r(pr_t$or_upper), p_holm = r(pr_t$p_holm), survives_holm = pr_t$survives_holm),
                  engaged = list(rate_current = r(pr_e$rate_if_1), rate_not = r(pr_e$rate_if_0), or = r(pr_e$odds_ratio), lo = r(pr_e$or_lower),
                                 hi = r(pr_e$or_upper), p = r(pr_e$p, 4), p_holm = r(pr_e$p_holm), survives_holm = pr_e$survives_holm)),
  frequency = list(took = list(rho = r(fq_t$rho), p_holm = r(fq_t$p_holm), p_fdr = r(fq_t$p_fdr), survives_fdr = fq_t$survives_fdr),
                   engaged = list(rho = r(fq_e$rho), p = r(fq_e$p, 4), p_holm = r(fq_e$p_holm))),
  sf8 = list(took_rho = r(bv("took", "sf8")$rho), engaged_rho = r(bv("engaged", "sf8")$rho)),
  model = list(predictors = unlist(pj$model$predictors), n = pj$model$complete_cases,
               practice_took = coef_of("took", "practice_current"), practice_engaged = coef_of("engaged", "practice_current"),
               sf8_took = coef_of("took", "z_sf8"), sf8_engaged = coef_of("engaged", "z_sf8"),
               auc_took = auc_of("took"), auc_engaged = auc_of("engaged")),
  main_wave = list(took = wave("took"), engaged = wave("engaged")),
  dose_rho_min = r(min(dose, na.rm = TRUE)), dose_rho_max = r(max(dose, na.rm = TRUE)),
  interactions_min_p_holm = r(min(inter$p_holm)),
  paused = list(n = pj$paused$n, on_page_median = r(dur$median, 0), on_page_q25 = r(dur$q25, 0), on_page_q75 = r(dur$q75, 0),
                ge60 = pj$paused$on_page_ge_60s, notes = pj$paused$wrote_note,
                sf8_rho = pcor("on-page seconds (paused)", "sf8")$rho,
                minutes_before = pcor("on-page seconds (paused)", "minutes_before")),
  declined = pj$decline_reasons$decliners,
  decline = lapply(Filter(function(z) !is.null(z[["pct", exact = TRUE]]), pj$decline_reasons$table),
                   function(z) list(reason = z$reason, n = z$n, pct = z[["pct", exact = TRUE]], sf8_mean = z$sf8_mean, tired_mean = z$tired_mean)),
  decline_sf8_p = r(pj$decline_reasons$sf8_by_reason_kruskal_all5),
  decline_tired_p = signif(pj$decline_reasons$tired_by_reason_kruskal_all5, 2)
)

# The 8 and FFMQ Observing (12) ---------
vs <- one(o("12_a_variance_share.csv"), direction == "8|obs" & label == "r2")
rc12 <- o("12_a_residual_correlations.csv")
kgl <- o("12_a_known_groups_latent.csv")
cp <- o("12_c_item_cross_pairs.csv")
bif <- one(o("12_b_bifactor_fit.csv"), grepl("reference", model))
bl <- o("12_b_bifactor_loadings.csv")[grepl("reference", model) & item %in% ld[form == 8, item]]
bl[, reversed := item %in% ld[form == 8 & content < 0, item]]
r12 <- function(s) { a <- one(rc12, score == s & grepl("^the 8", focal)); b <- one(rc12, score == s & grepl("^FFMQ", focal))
  list(score = s, label = lab(s), total = r(a$r_total), net8 = r(a$residual_r), net8_lo = r(a$lower), net8_hi = r(a$upper), net_obs = r(b$residual_r)) }
overlap <- list(
  r2 = r(vs$est), lo = r(vs$lower), hi = r(vs$upper),
  cross_max = { x <- cp[which.max(abs(keyed_residual))]; list(statement = x$item_8, code = one(ld, form == 8 & item == x$item_8)$source,
                                                             residual = r(x$residual), keyed = r(x$keyed_residual)) },
  rows = lapply(c("sbi_moment", "bfi_aesthetic", "maia_noticing", "pss4", "phq4_total", "lifesat", "crit_pause_days", "ncs6", "bidr_sde"), r12),
  practice = list(total = r(one(kgl, model == "raw" & label == "d_8")$est),
                  net8 = r(one(kgl, model == "8|obs" & label == "d_residual")$est), net8_lo = r(one(kgl, model == "8|obs" & label == "d_residual")$lower),
                  net8_hi = r(one(kgl, model == "8|obs" & label == "d_residual")$upper),
                  net_obs = r(one(kgl, model == "obs|8" & label == "d_residual")$est)),
  bifactor = list(specific_share = r(bif$specific_share_of_reliable_8), omega_s = r(bif$omega_s_specific_8), cfi = r(bif$cfi.scaled),
                  s_reversed = r(range(abs(bl[reversed == TRUE, S8]))), s_forward = r(range(abs(bl[reversed == FALSE, S8]))))
)

# All four H2 neighbours (18) ---------
rg <- o("18_regression.csv"); rr <- o("18_residual_reliability.csv"); rc18 <- o("18_residual_correlations.csv"); p18 <- o("18_practice.csv")
q18 <- function(v, q) one(rg, variant == v & quantity == q)
pq <- function(v, a, q) one(p18, variant == v & analysis == a & quantity == q)
r18 <- function(s) { x <- one(rc18, score == s)
  list(score = s, label = lab(s), total = r(x$r_total), resid = r(x$residual_r), lo = r(x$lower), hi = r(x$upper), resid_af = r(x$residual_r_af)) }
residual4 <- list(
  r2 = r(q18("registered", "r2")$est), lo = r(q18("registered", "r2")$lower), hi = r(q18("registered", "r2")$upper),
  r2_af = r(q18("agreement_free", "r2")$est),
  beta = lapply(c("obs", "sbi", "aes", "not"), function(k) { x <- q18("registered", paste0("beta_", k)); list(label = lab(x$neighbour), beta = r(x$est), lo = r(x$lower), hi = r(x$upper)) }),
  reliability = r(one(rr, variant == "registered" & grepl("omega", reliability_of_8))$residual_reliability),
  reliability_af = r(one(rr, variant == "agreement_free" & grepl("omega", reliability_of_8))$residual_reliability),
  rows = lapply(c("vig_self", "crit_pause_days", "maia_bodylistening", "maia_selfregulation", "spsq_sd", "ffmq_acting", "pss4",
                  "ncs6", "bidr_sde", "bidr_im", "lifesat", "engaged_pause", "separation"), r18),
  practice_group = list(total = r(pq("registered", "practice group", "d_8_whole")$est), resid = r(pq("registered", "practice group", "d_residual")$est),
                        lo = r(pq("registered", "practice group", "d_residual")$lower), hi = r(pq("registered", "practice group", "d_residual")$upper)),
  practice_dose = list(total = r(pq("registered", "practice frequency (dose)", "slope_8_whole")$est),
                       resid = r(pq("registered", "practice frequency (dose)", "slope_residual")$est),
                       lo = r(pq("registered", "practice frequency (dose)", "slope_residual")$lower),
                       hi = r(pq("registered", "practice frequency (dose)", "slope_residual")$upper)),
  observing_af = r(q18("agreement_free", "totnb_obs")$est)
)

# A general factor across the battery (13) ---------
gj <- oj("13_exploratory_general_factor.json"); gl <- o("13_general_loadings.csv"); ng <- o("13_net_of_general.csv")
am <- matrix(unlist(gj$agreement_between_methods), nrow = length(gj$agreement_between_methods), byrow = TRUE)
g13 <- function(s) { x <- one(ng, g_definition == "sl" & measure == s); list(score = s, label = lab(s), r = r(x$r_with_8), share = r(x$share_accounted), net = r(x$partial_r_net_of_G)) }
pgn <- gj$practice_group_net
general <- list(
  k = gj$schmid_leiman$k, measures = length(gj$matrix$battery), share_common = r(gj$schmid_leiman$g_variance_share_of_common),
  alternatives = nrow(am) - 1L, agreement_min = r(min(am)), smoothed = gj$matrix$smoothed, heywood = lab(gj$heywood$primary),
  loadings = rows(gl[, .(score = measure, label = lab(measure), g = r(sl_g))]),
  net = lapply(c("bidr_sde", "eq_decentering", "ncs6", "ffmq_observing", "crit_pause_days"), g13),
  wellbeing_net = lapply(c("pss4", "lifesat", "phq4_total"), g13),
  n_ge_40 = gl[abs(sl_g) >= .40, .N],
  heywood_max_shift = r(max(abs(gl$sl_g - gl$sl_g_without_bfi_intellectual), na.rm = TRUE)),
  practice = list(d_latent = r(pgn$d_07_latent), net_min = r(min(vapply(pgn$net, `[[`, 0, "d_net_of_G"))), net_max = r(max(vapply(pgn$net, `[[`, 0, "d_net_of_G"))))
)

# Statements (14) and swapping U4 (15) ---------
loo <- o("14_leave_one_out.csv")[dropped != "none"]
ij <- oj("14_exploratory_items.json")
spec <- rbindlist(lapply(ij$net_of_content, function(z) data.table(af = z$agreement_free, target = z$target, code = z$source, role = z$role, r = z$specific_r_keyed)))
spec8 <- spec[target == "ffmq_observing" & role == "the 8" & af == FALSE]
rev8 <- ld[form == 8 & content < 0, source]
afv <- one(o("04_agreement_factor_validity.csv"), form == "sf8" & score == "ffmq_observing")
sb <- o("15_side_by_side.csv")
sbq <- function(q) as.list(one(sb, quantity == q)[, .(registered, swap, seven)])
items <- list(
  loo = list(obs_min = r(min(loo$r_ffmq_observing)), obs_max = r(max(loo$r_ffmq_observing)), ncs_min = r(min(loo$r_ncs6)), ncs_max = r(max(loo$r_ncs6)),
             sde_min = r(min(loo$r_bidr_sde)), sde_max = r(max(loo$r_bidr_sde)), all_meet_h1 = all(loo$h1_rule_met)),
  specific_reversed = r(range(spec8[code %in% rev8, r])), specific_forward = r(range(spec8[!code %in% rev8, r])),
  specific_af_max = r(max(abs(spec[target == "ffmq_observing" & role == "the 8" & af == TRUE, r]))),
  agreement_with_observing = r(afv$r_with_agreement),
  observing_af = residual4$observing_af
)
f15 <- o("15_fit.csv")   # omega and alpha at 4 decimals (15_side_by_side.csv rounds them to 3)
f15q <- function(col) setNames(lapply(c("registered", "swap", "seven"), function(k) r(f15[set == k][[col]])), c("registered", "swap", "seven"))
swap <- list(cfi = sbq("CFI"), rmsea = sbq("RMSEA"), h1 = sbq("H1 rule met"), omega = f15q("omega"), alpha = f15q("alpha"),
             r_obs = sbq("latent r: ffmq_observing"), r_ncs = sbq("latent r: ncs6"), r_sde = sbq("latent r: bidr_sde"), d = sbq("d: H6: current vs never"))
swap[c("cfi", "rmsea", "h1", "r_obs", "r_ncs", "r_sde", "d")] <- lapply(swap[c("cfi", "rmsea", "h1", "r_obs", "r_ncs", "r_sde", "d")], function(x) lapply(x, function(v) {
  v <- sub(" \\[.*", "", as.character(v)); if (v %in% c("TRUE", "FALSE")) as.logical(v) else as.numeric(v) }))

# Practice as a dose (16) ---------
dl <- o("16_dose_levels.csv")[code == "freq5" & outcome == "the 8 (observed)"]
dp <- o("16_pause.csv")[by == "freq5"]
tr <- o("16_trend.csv"); gr <- o("16_groups.csv"); sp <- o("16_spirituality.csv"); ty <- o("16_types.csv")
gp <- o("16_pause.csv")[by == "group"]
practice <- list(
  levels = lapply(seq_len(nrow(dl)), function(i) list(level = dl$level[i], label = sub("^[0-9]: ", "", dl$label[i]), n = dl$n[i], mean = r(dl$mean[i]),
                                                      sd = r(dl$sd[i]), took_pct = dp[level == dl$level[i], took_pct])),
  slope = r(one(tr, analysis == "the 8 on freq5" & sample == "all")$slope), rho = r(one(tr, analysis == "the 8 on freq5" & sample == "all")$spearman_rho),
  nonlinear_p = r(one(tr, analysis == "the 8 on freq5" & sample == "all")$nonlinearity_p),
  slope_spirit = r(one(tr, analysis == "the 8 on freq5, spirituality and religiosity held constant" & sample == "all")$slope),
  spiritual_rho = r(one(sp, variable == "spiritual")$spearman_rho_with_8), religious_rho = r(one(sp, variable == "religious")$spearman_rho_with_8),
  net_g_obs = { x <- one(tr, analysis == "the 8 on freq5, net of G_obs and Observing"); list(slope = r(x$slope), lo = r(x$lower), hi = r(x$upper)) },
  dedicated_vs_current = list(d = r(one(gr, contrast == "dedicated vs current, not dedicated" & covariates == "")$d),
                              d_spirit = r(one(gr, contrast == "dedicated vs current, not dedicated" & covariates != "")$d)),
  yoga = list(diff = r(one(ty, type == "yoga")$diff_adj_freq), d = r(one(ty, type == "yoga")$d_adj_freq)),
  dedicated_pause = list(took = one(gp, level == "dedicated")$took, n = one(gp, level == "dedicated")$n)
)

# The vignettes (17) ---------
vj <- oj("17_exploratory_vignettes.json"); vi <- o("17_vignette_items.csv")[kind != "self-rating"]; mo <- o("17_moderation.csv"); sr <- o("17_self_rating.csv")
amb <- vi$vignette %in% unlist(vj$flagged_ambiguous)
vignettes <- list(
  ambiguous = sub("vig_", "", unlist(vj$flagged_ambiguous)),
  r_rest_ambiguous = r(range(vi[amb, r_rest_keyed])), r_rest_others = r(range(vi[!amb, r_rest_keyed])),
  without = list(tells_apart = vj$rescore$without_2_4_7$tells_apart, share = r(vj$rescore$without_2_4_7$share), n = vj$rescore$without_2_4_7$n),
  max_abs_r_8_or_self = r(max(abs(c(vi$r_with_8, vi$r_with_self)))),
  self_mean = list(no = r(one(sr, tells_apart == "no")$self_mean), yes = r(one(sr, tells_apart == "yes")$self_mean)),
  contrasts = nrow(mo), holding_holm = sum(mo$wald_p_holm < .05),
  moderation = rows(mo[target %in% c("ncs6", "bidr_sde", "crit_pause_days", "ffmq_observing"),
                       .(target, label = lab(target), no = r(est_no), yes = r(est_yes), diff = r(difference_yes_minus_no), lo = r(diff_lower), hi = r(diff_upper),
                         p_holm = r(wald_p_holm))])
)

# Agreement across the battery (20) ---------
# Read from 20's JSON, not its CSVs: the CSVs round to 3 decimals, which the page would round again (.7045 shown as .71).
# The registered and ri_free columns carry 04/12/13/18's own 3-decimal values either way.
aj <- fromJSON(sf2b_out("20_exploratory_battery_agreement.json"), simplifyVector = TRUE)
tb <- function(x) as.data.table(x)
ac <- tb(aj$a$correlations); ah <- tb(aj$a$h3b); ab <- tb(aj$b$h3d)
avs <- tb(aj$c$variance_share); arc <- tb(aj$c$residual_correlations)
arg <- tb(aj$d$regression); ard <- tb(aj$d$residual_correlations); apr <- tb(aj$d$practice_group)
aes <- tb(aj$e$summary); agl <- tb(aj$e$loadings); ang <- tb(aj$e$net_of_G); apg <- tb(aj$e$practice_group_net)
afd <- tb(aj$f$free_correlation_diagnostic); asz <- tb(aj$f$size); arv <- tb(aj$reversed_items_by_score)
a20 <- function(s) { x <- one(ac, score == s)
  list(score = s, label = lab(s), registered = r(x$r_registered), ri_free = r(x$r_ri_free), battery = r(x$r_battery_A), lo = r(x$lower), hi = r(x$upper)) }
c20 <- function(s) { x <- one(arc, score == s)
  list(score = s, label = lab(s), registered = r(x$residual_r_registered), battery = r(x$residual_r_battery_A), lo = r(x$lower), hi = r(x$upper)) }
d20 <- function(s) { x <- one(ard, score == s)
  list(score = s, label = lab(s), registered = r(x$residual_r_registered), agreement_free = r(x$residual_r_agreement_free_18),
       battery = r(x$residual_r_battery_A), lo = r(x$lower), hi = r(x$upper)) }
h3d20 <- function(q) one(ab, quantity == q)
es20 <- function(m) one(aes, matrix_source == m)
ng20 <- function(m, s) r(one(ang, matrix_source == m & g_definition == "sl" & measure == s)$partial_r_net_of_G)
pg20 <- function(m, g = "sl") r(one(apg, matrix_source == m & g_definition == g)$d_net_of_G)
vs20 <- one(avs, quantity == "unexplained"); un4 <- one(arg, quantity == "unexplained"); dres <- one(apr, quantity == "d_residual")
free_sig <- afd[usable == TRUE & p < .05][order(-abs(r_A_with_factor))]
agreement <- list(
  n = aj$n,
  correlations = lapply(c("ffmq_observing", "maia_noticing", "sbi_moment", "bfi_aesthetic", "samd_md", "eq_decentering",
                          "bidr_sde", "ncs6", "bsss8", "asi_pc"), a20),
  h3b = lapply(seq_len(nrow(ah)), function(i) list(contrast = ah$contrast[i], first = lab(ah$first[i]), second = lab(ah$second[i]),
                                                     registered = r(ah$difference_registered_refit[i]), battery = r(ah$difference_battery_A[i]),
                                                     lo = r(ah$lower[i]), hi = r(ah$upper[i]), above_zero = ah$above_zero[i])),
  h3d = list(registered = r(h3d20("R registered (04, bootstrap CI)")$est), registered_lo = r(h3d20("R registered (04, bootstrap CI)")$lower),
             registered_hi = r(h3d20("R registered (04, bootstrap CI)")$upper),
             battery = r(h3d20("R battery A (bootstrap CI, 04's resamples)")$est), lo = r(h3d20("R battery A (bootstrap CI, 04's resamples)")$lower),
             hi = r(h3d20("R battery A (bootstrap CI, 04's resamples)")$upper),
             resamples_ok = h3d20("R battery A (bootstrap CI, 04's resamples)")$resamples_ok,
             resamples = h3d20("R battery A (bootstrap CI, 04's resamples)")$resamples,
             upper_below_80 = h3d20("R battery A (bootstrap CI, 04's resamples)")$upper_below_80,
             agreement_free = r(h3d20("R agreement free (18's agreement_free form, delta CI)")$est)),
  after_observing = list(unexplained_registered = r(vs20$est_registered), unexplained = r(vs20$est_battery_A), lo = r(vs20$lower), hi = r(vs20$upper),
                         rows = lapply(c("sbi_moment", "bfi_aesthetic", "pss4", "crit_pause_days", "ncs6", "bidr_sde"), c20)),
  after_four = list(unexplained_registered = r(un4$est_registered), unexplained_agreement_free = r(un4$est_agreement_free_18),
                    unexplained = r(un4$est_battery_A), lo = r(un4$lower), hi = r(un4$upper),
                    rows = lapply(c("vig_self", "crit_pause_days", "spsq_sd", "maia_bodylistening", "ffmq_nonreactivity", "ncs6", "bidr_sde"), d20),
                    practice = list(n = dres$n, registered = r(dres$est_registered), agreement_free = r(dres$est_agreement_free_18),
                                    battery = r(dres$est_battery_A), lo = r(dres$lower), hi = r(dres$upper), p = r(dres$p))),
  general = list(group_factors = es20("battery A")$group_factors, smallest_eigenvalue = r(es20("battery A")$smallest_eigenvalue),
                 smallest_eigenvalue_13 = r(es20("13 (registered agreement, reproduction check)")$smallest_eigenvalue),
                 heywood = lab(strsplit(es20("battery A")$heywood, " ")[[1]]),
                 share_common = r(es20("battery A")$g_share_of_common),
                 g8_registered = r(one(agl, measure == "sf8")$sl_g_13), g8 = r(one(agl, measure == "sf8")$sl_g_battery_A),
                 g8_calibrated = r(one(agl, measure == "sf8")$sl_g_calibrated),
                 g_obs_registered = r(one(agl, measure == "ffmq_observing")$sl_g_13), g_obs = r(one(agl, measure == "ffmq_observing")$sl_g_battery_A),
                 obs_net_registered = ng20("13", "ffmq_observing"), obs_net = ng20("battery A", "ffmq_observing"),
                 obs_net_calibrated = ng20("battery A, calibrated", "ffmq_observing"),
                 practice_net_registered = pg20("13"), practice_net = pg20("battery A"), practice_net_calibrated = pg20("battery A, calibrated"),
                 practice_net_min = r(min(apg[matrix_source == "battery A", d_net_of_G])), practice_net_max = r(max(apg[matrix_source == "battery A", d_net_of_G]))),
  agreement_share_8 = r(one(asz, instrument == "sf8")$a_share_mean),
  free = lapply(seq_len(nrow(free_sig)), function(i) list(score = free_sig$score[i], label = lab(free_sig$score[i]),
                                                            r = r(free_sig$r_A_with_factor[i]), lo = r(free_sig$lower[i]), hi = r(free_sig$upper[i]))),
  free_maia_unusable = all(!afd[grepl("^maia_", score), usable]),   # the Brief MAIA-2 diagnostic gave no usable fit
  all_reversed = lab(arv[reversed == items & reversed > 0, score]),
  none_reversed_neighbours = lab(arv[reversed == 0 & score %in% c("ffmq_observing", "maia_noticing"), score])
)

# Subsets of the 13 statements (21) ---------
sj <- oj("21_exploratory_subsets.json")
bs <- o("21_best_subsets.csv"); bt <- o("21_bootstrap_test_half.csv"); ss <- o("21_stability_summary.csv")
sst <- o("21_stability_splits.csv"); ssub <- o("21_stability_subsets.csv"); sstm <- o("21_stability_statements.csv"); sch <- o("21_stability_chosen.csv")
bench <- function(role_name) {
  f <- function(v) one(bs, role == role_name & sample == "full" & version == v)
  x <- f("i"); y <- f("ii")
  list(statements = x$statements, k = x$k, omega = r(x$omega), r_obs = r(x$r_obs), r_obs_af = r(y$r_obs), R = r(x$R_h3d),
       pause = r(x$r_pause_days), vig = r(x$r_vig_self), h6_d = r(x$h6_d))
}
p1 <- sst[split == 1]
chosen <- lapply(seq_len(nrow(p1)), function(i) {
  x <- p1[i]; b <- one(bt, version == x$version & objective == x$objective)
  t <- one(bs, id == b$id & sample == "test" & version == x$version)
  list(version = x$version, objective = x$objective, statements = x$statements, k = x$k,
       search = r(x$chosen_search), test = r(x$chosen_test), the8_search = r(x$the8_search), the8_test = r(x$the8_test),
       gain_search = r(x$gain_search), gain_test = r(x$gain_test), boot_lo = r(b$difference_lower), boot_hi = r(b$difference_upper),
       boot_used = b$used, boot_resamples = b$resamples, test_h1 = x$test_h1_met, test_omega = r(x$test_omega), test_constraints_met = x$test_proper,
       convergent_flag = t$convergent_flag, test_h6_d = r(t$h6_d), material = x$material_in_test)
})
top5 <- function(v, objs) {   # the five best per objective in the primary split, full-sample rows
  pat <- paste0("best [1-5] on (", paste(objs, collapse = "|"), "), version ", v, "($|;)")
  bs[sample == "full" & version == v & grepl(pat, role)]
}
reg_top <- top5("i", c("r_obs", "R")); af_ch <- sch[version == "ii" & rank == 1]
share_of <- function(v, obj, code) r(one(sstm, version == v & objective == obj & statement == code & top == 1)$share)
# The page shows the nearest miss at 3 decimals, so it is rounded here at 3 (4 then 3 would round -.04448 to -.045).
miss <- function(v) { x <- one(bs, statements == "N4 sfc_36 sfc_35 N6 sf6_44" & sample == "test" & version == v)
  list(gain = r(x$d_r_obs, 3), omega = r(x$omega), h1 = x$h1_rule_met, flag = x$convergent_flag, proper = !isTRUE(x$improper_validity)) }
subsets <- list(
  n = sj$n, candidates = sj$candidates, splits = sj$splits, half_n = sj$primary_split$search_n,
  material_primary = length(sj$material_in_primary_test_half),
  benchmarks = list(the8 = bench("benchmark: the_8"), the10 = bench("benchmark: the_10"), the7 = bench("benchmark: the_7_without_U4")),
  chosen = chosen,
  stability = lapply(seq_len(nrow(ss)), function(i) { x <- ss[i]
    list(version = x$version, objective = x$objective, distinct = x$distinct_subsets, max_chosen = max(ssub[version == x$version & objective == x$objective, times_chosen]),
         gain_search = r(x$mean_gain_search), gain_test = r(x$mean_gain_test), test_omega_ok = x$test_omega_ok) }),
  max_chosen = max(ssub$times_chosen),
  registered_winners = list(h6_min = r(min(reg_top$h6_d)), h6_max = r(max(reg_top$h6_d)), pause_min = r(min(reg_top$r_pause_days)),
                            pause_max = r(max(reg_top$r_pause_days)),
                            share = lapply(c("sf6_36", "sfc_x2"), function(s) list(code = s, r_obs = share_of("i", "r_obs", s), R = share_of("i", "R", s)))),
  af_winners = list(k_min = min(af_ch$k), k_max = max(af_ch$k), omega_median = r(median(af_ch$omega)),
                    # 21's test_proper is every constraint: the H1 rule, omega >= .75, proper, oriented as keyed
                    test_constraints_failed = sum(!sst[version == "ii", test_proper]), test_n = nrow(sst[version == "ii"]),
                    share = lapply(c("N4", "sfc_35", "sf6_44", "sfc_36"), function(s) list(code = s, r_obs = share_of("ii", "r_obs", s), R = share_of("ii", "R", s)))),
  nearest_miss = list(statements = "N4 sfc_36 sfc_35 N6 sf6_44", k = 5L, registered = miss("i"), agreement_free = miss("ii"))
)

# Write ---------
out <- list(
  generated = as.character(Sys.Date()),
  commit_note = "Analysis code: Normega/sense-foraging-scale, study2b/analysis (registered 01-09; deviation 10; exploratory 11-18, 20 and 21; this file 19).",
  recorded = recorded, sample = sample, decisions = decisions, h1 = h1, validity = validity, h4 = h4, h5 = h5, h6 = h6, h8 = h8,
  pause = pause, overlap = overlap, residual4 = residual4, general = general, items = items, swap = swap, practice = practice, vignettes = vignettes,
  agreement = agreement, subsets = subsets
)
writeLines(toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null", null = "null"), dest, useBytes = TRUE)
sf2b_log("wrote ", dest)
