# 08: the exploratory analyses (Section 5, "Exploratory analysis"). None is a decision.
#   H8       invariance of the H1 model across people who do and do not tell the modes
#            apart: configural, thresholds, loadings (Wu & Estabrook 2016), agreement
#            loadings fixed at 1 in both groups and its variance free; a step is
#            non-invariant when delta CFI < -.010 and delta RMSEA > .015
#   also     the same by gender (men, women) and by country group (if each has 200+)
#   trial    the 8 against the trial's Access + Valuing score in its own wording
#   pause    seconds on the page, log(1 + x), among those who paused; acceptance alone;
#            the 8 in decliners giving "not interested" against "no time" or "not a
#            good place"; take-up by tiredness, setting and device; the 8's association
#            with the engaged pause beside tiredness's
#   order    a measure's position within its half against its mean and its correlation with the 8
#   vignettes  separation, the share who tell the modes apart, the self-rating, each vignette
#   BIDR     dichotomous scoring for H3c (observed correlations)
# PREREG: Section 5 lists H8's steps as "configural, then loadings, then thresholds";
# for ordinal items Wu & Estabrook (2016) require thresholds before loadings, which is
# the order used here (and in Study 2's code). The text should be corrected before filing.

source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")
A <- readRDS(sf2b_derived("sf2b")); D <- A$defs; M <- A$meta; S <- A$scores; X <- A$scale_raw
out <- list()

# ---- invariance --------------------------------------------------------------------------------
inv <- list(sf2b_run_mg(X, D$sf8, M$tells_apart, "H8: tells the modes apart (yes / no)"),
            sf2b_run_mg(X, D$sf8, M$gender2, "gender (man / woman)"))
if (isTRUE(A$country_rule$analysed)) inv[[3]] <- sf2b_run_mg(X, D$sf8, M$country_group, "country group")
# A sequence that stopped at a configural model that failed has no step to judge: give its
# table the column, empty, so the summaries below still run.
inv <- lapply(inv, function(r) { if (!"noninvariant" %in% names(r$table)) r$table[, noninvariant := NA]; r })
fwrite(rbindlist(lapply(inv, `[[`, "table"), fill = TRUE), sf2b_out("08_invariance.csv"))
mg <- rbindlist(lapply(inv, `[[`, "merges")); if (nrow(mg)) fwrite(mg, sf2b_out("08_invariance_category_merges.csv"))
out$invariance <- lapply(inv, function(r) list(groups = r$groups, released = r$free,
                                               last_step = tail(r$table$step, 1), noninvariant_steps = r$table[noninvariant == TRUE, step]))
out$country_rule <- A$country_rule

# ---- the trial's measure ----------------------------------------------------------------------------
out$trial <- list(sf8_with_trial_av10 = sf2b_r_ci(S$sf8, S$trial_av10), sf10_with_trial_av10 = sf2b_r_ci(S$sf10, S$trial_av10),
                  note = "each one's correlations with the validity measures: 04_observed_validity.csv (score_of = trial_av10)")

# ---- the pause, further -----------------------------------------------------------------------------
P <- M[!is.na(engaged_pause)]; SP <- S[match(P$participant, S$participant)]
took <- P$pause_taken == 1
reason <- P$pause_decline_reason
y_dec <- fifelse(reason %in% "not_interested", 1L, fifelse(reason %in% c("no_time", "not_good_place"), 0L, NA_integer_))
z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
tired <- P$pause_tired
lg <- glm(engaged_pause ~ zs + zt, family = binomial, data = data.frame(engaged_pause = P$engaged_pause, zs = z(SP$sf8), zt = z(tired)))
rate_by <- function(v) { d <- data.table(level = as.character(v), taken = P$pause_taken, engaged = P$engaged_pause)
  d[, .(n = .N, took = r3(mean(taken)), engaged = r3(mean(engaged))), by = level][order(level)] }
out$pause <- list(
  seconds_on_page_log1p_with_sf8 = sf2b_r_ci(SP$sf8[took], log1p(P$pause_on_page_s[took]), "spearman"),
  accepted_with_sf8 = sf2b_r_ci(SP$sf8, P$pause_taken, "spearman"),
  not_interested_vs_no_time_or_place = if (length(unique(na.omit(y_dec))) == 2) sf2b_kg(SP$sf8, y_dec) else NULL,
  engaged_on_sf8_and_tiredness = as.data.table(summary(lg)$coefficients, keep.rownames = "term"),
  by_tiredness = rate_by(tired), by_setting = rate_by(P$pause_setting), by_device = rate_by(P$pause_device),
  decline_reasons = as.list(table(reason[!took], useNA = "ifany")),
  ended_by = as.list(table(P$pause_ended_by[took], useNA = "ifany")),
  notes_written = sum(P$pause_has_note, na.rm = TRUE), notes = "coded blind to scores; the text is in derived/, never in output/")

# ---- order effects -------------------------------------------------------------------------------------
pos <- A$validity_positions
slug_scores <- c(lapply(D$validity, function(v) setNames(names(v$scores), rep(v$slug, length(v$scores)))))
pairs <- rbindlist(lapply(D$validity, function(v) data.table(slug = v$slug, score = names(v$scores))))
pairs <- rbind(pairs, data.table(slug = "lifesat", score = "lifesat"),
               data.table(slug = "sf2-past-week", score = c("crit_pause_days", "crit_notice_first", "crit_outdoors")))
ord <- rbindlist(lapply(seq_len(nrow(pairs)), function(k) {
  sl <- pairs$slug[k]; sc <- pairs$score[k]
  if (!(sl %in% names(pos))) return(NULL)
  p <- pos[[sl]][match(S$participant, pos$participant)]; y <- S[[sc]]
  early <- !is.na(p) & p <= median(p, na.rm = TRUE)
  data.table(measure = sl, score = sc, n = sum(!is.na(p) & !is.na(y)),
             rho_position_score = r3(sf2b_r_ci(p, y, "spearman")$r), p = signif(sf2b_r_ci(p, y, "spearman")$p, 3),
             r_with_sf8_early = r3(sf2b_r_ci(S$sf8[early], y[early])$r), r_with_sf8_late = r3(sf2b_r_ci(S$sf8[!early & !is.na(p)], y[!early & !is.na(p)])$r))
}))
fwrite(ord, sf2b_out("08_order_effects.csv"))

# ---- vignettes ----------------------------------------------------------------------------------------
ta <- M$tells_apart; k <- sum(ta == "yes", na.rm = TRUE); nn <- sum(!is.na(ta))
wil <- if (nn) prop.test(k, nn, correct = FALSE)$conf.int else c(NA, NA)
out$vignettes <- list(separation_mean = r3(mean(M$separation, na.rm = TRUE)), separation_sd = r3(sd(M$separation, na.rm = TRUE)),
                      n_with_separation = nn, tells_apart = k, share = r3(k / nn), share_ci_wilson = r3(wil),
                      self_rating_with_sf8 = sf2b_r_ci(S$sf8, M$vig_self),
                      vignette_means = as.list(round(colMeans(A$vignettes[, names(D$vkind)], na.rm = TRUE), 2)),
                      vignette_kind = as.list(D$vkind))

# ---- BIDR dichotomous ----------------------------------------------------------------------------------
out$bidr_dichotomous <- list(sde = sf2b_r_ci(S$sf8, S$bidr_sde_dichotomous), im = sf2b_r_ci(S$sf8, S$bidr_im_dichotomous))

rnd <- function(x) if (is.data.frame(x)) x else if (is.list(x)) lapply(x, rnd) else if (is.numeric(x)) signif(x, 4) else x
write_json(rnd(out), sf2b_out("08_exploratory.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(rbindlist(lapply(inv, `[[`, "table"), fill = TRUE)[, .(analysis, step, released, cfi, rmsea, d_cfi, d_rmsea, noninvariant)])
