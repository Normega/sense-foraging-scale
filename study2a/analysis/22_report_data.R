# 22: the data behind senseforaging.com/research/study2a/results: aggregates only, read from
# the outputs of 01-21 (real mode). Writes study2aResults.json to SF2A_REPORT_OUT (default:
# the Handbook checkout's src/pages/research/). Rerun after any of 01-21 changes.

source("study2a/analysis/config.R")
stopifnot(SF2A$mode == "real")
o <- function(f) fread(sf2a_out(f))
r <- function(x, d = 3) round(x, d)
def <- sf2a_definition()
dest <- Sys.getenv("SF2A_REPORT_OUT", "F:/gits/Handbook/src/pages/research/study2aResults.json")

flow <- o("01_flow.csv"); flags <- o("01_flags.csv")[excluded_attention == FALSE]
vd <- o("03_vignettes_decision.csv"); ve <- o("03_vignettes_each.csv")
vj <- fromJSON(SF2A$items, simplifyVector = FALSE)
vt <- setNames(vapply(vj$vignettes, function(v) v$kind, ""), sprintf("vig_%s", sub("^V", "", vapply(vj$vignettes, `[[`, "", "id"))))
vx <- setNames(vapply(vj$vignettes, function(v) v$text, ""), names(vt))
pa <- o("04_parallel_analysis.csv"); fac <- o("05_factors_k04.csv"); rel <- o("08_reliability_k04.csv")
fit5 <- o("05_primary_fit_k04.csv"); ri <- o("05_agreement_factor_k04.csv"); ius <- o("07_ius.csv"); L <- o("06_decision_L_k04.csv")
fc <- as.matrix(o("05_factor_correlations_k04.csv")[, -1])
sens <- o("10_sensitivity_k04.csv")
pe <- o("12_plain_efa_factors.csv"); ev <- o("12_eigenvalues.csv"); gfit <- o("12_method_models_fit.csv"); gri <- o("12_G_agreement.csv")
av <- jsonlite::fromJSON("F:/radlab_data/sf_study1/output/real/12_access_valuing_study1.json")
tri <- o("14_chosen_summary.csv"); al <- o("15_alpha.csv"); st <- o("15_structure.csv"); cg <- o("15_efa_congruence.csv"); gm <- o("15_G_model.csv")
bg <- o("16_by_background.csv"); sr <- o("16_separation_regression.csv"); sc <- o("16_scores_by_third.csv"); scc <- o("16_score_correlations.csv")
forms <- o("21_forms.csv"); cfa <- o("21_cfa.csv"); stm <- o("21_statements.csv"); inv <- o("21_invariance_preview.csv")
mf <- o("19_mode_free.csv"); c20 <- o("20_fit.csv")
labels <- c(F1 = "Pull of the doing mode", F2 = "The modes are the same", F3 = "Broad sensing", F4 = "Need to explain or act")

out <- list(
  generated = as.character(Sys.Date()),
  commit_note = "Analysis code: Normega/sense-foraging-scale, study2a/analysis (registered 01-11; exploratory 12-21).",
  sample = list(enrolled = flow[step == "enrolled (not staff tests)", n], submitted = flow[step == "submitted the scale", n],
                excluded = flow[step == "excluded: failed two or more attention checks", n], analysable = flow[step == "analysable", n],
                flagged = flags$any_flag, one_check = flags$one_check, speeder = flags$speeder, longstring = flags$longstring,
                practice = setNames(as.list(o("11_practice_groups.csv")[!is.na(practice_group) & practice_group != "", N]),
                                    o("11_practice_groups.csv")[!is.na(practice_group) & practice_group != "", practice_group])),
  vignettes = list(pct = vd$pct_tell_apart, lo = vd$pct_ci_lo, hi = vd$pct_ci_hi, n_tell = vd$n_tell_apart, threshold = vd$threshold_pct,
                   sensing = vd$mean_sensing, doing = vd$mean_doing, midpoint = vd$midpoint, separation = vd$mean_separation,
                   each = lapply(seq_len(nrow(ve)), function(i) list(id = sub("vig_", "V", ve$vignette[i]), title = unname(vt[ve$vignette[i]]),
                                                                    text = unname(vx[ve$vignette[i]]), kind = ve$kind[i],
                                                                    mean = r(ve$mean[i], 2), wrong_side = ve$wrong_side[i]))),
  structure = list(pa = pa[, which(real < random_95)[1] - 1L], eigen = r(ev$eigenvalue[1:7], 1), first_pct = ev$pct[1],
                   fit = as.list(fit5[1, .(cfi = r(cfi.scaled), rmsea = r(rmsea.scaled), srmr = r(srmr))]), agreement_pct = ri$mean_pct_item_variance,
                   factors = lapply(seq_len(nrow(fac)), function(i) list(f = fac$factor[i], label = unname(labels[fac$factor[i]]), n = fac$n_primary_ge_40[i],
                                                                         n_rev = fac$n_reversed_ge_40[i], composition = fac$composition[i],
                                                                         omega = r(rel$ordinal_omega[i], 2))),
                   r_F1_F3 = r(fc[1, 3], 2), r_F3_F4 = r(fc[3, 4], 2),
                   letting_be_forward_ge40 = L$n_forward_ge_40, letting_be_forward = L$n_forward,
                   ius = list(forward = r(ius$r_ius_forward, 2), reversed = r(ius$r_ius_reversed, 2), t2 = r(ius$t2, 2)),
                   congruence = list(tell_apart = range(sens[analysis == "tell_apart", congruence]),
                                     without_flagged = range(sens[analysis == "without_flagged", congruence]))),
  followup = list(
    plain_efa = lapply(1:5, function(k) lapply(which(pe$k == k), function(i) list(n = pe$n_ge_40[i], n_rev = pe$n_reversed[i], composition = pe$composition[i]))),
    G = list(cfi = r(gfit[model == "G", cfi.scaled]), rmsea = r(gfit[model == "G", rmsea.scaled]), agreement_pct = gri[factor == "RI, G model", mean_pct_item_variance]),
    access_valuing = list(s1_r = av$AV10$two_factor$latent_r, s2a_r = .793, s2a_lo = .738, s2a_hi = .847,
                          s1_cfi2 = av$AV10$two_factor$fit$cfi.scaled, s1_cfi1 = av$AV10$one_factor$fit$cfi.scaled, s2a_cfi2 = .943, s2a_cfi1 = .911),
    triples = lapply(seq_len(nrow(tri)), function(i) list(factor = tri$factor[i], statements = tri$statements[i], omega_sub = r(tri$omega_sub[i], 2),
                                                          specific = r(tri$omega_specific[i], 2), specific_validated = r(tri$omega_specific_validated[i], 2),
                                                          alpha = r(tri$alpha[i], 2))),
    vignette_groups = list(alpha = lapply(c("best 8 (12)", "trial's 10, Study 1 wording", "all 58"), function(s) {
                             a <- al[contrast == "registered" & set == s]; b <- al[contrast == "thirds" & set == s]
                             list(set = s, hi = r(a$alpha_hi, 2), lo = r(a$alpha_lo, 2), top = r(b$alpha_hi, 2), bottom = r(b$alpha_lo, 2)) }),
                           r_fwd_rev = list(hi = r(st[contrast == "registered" & group == "hi", r_composites], 2), lo = r(st[contrast == "registered" & group == "lo", r_composites], 2)),
                           agreement = list(hi = r(gm[contrast == "registered" & group == "hi", agreement_pct_item_variance], 1),
                                            lo = r(gm[contrast == "registered" & group == "lo", agreement_pct_item_variance], 1)),
                           G_rev = list(hi = r(gm[contrast == "registered" & group == "hi", G_mean_keyed_reversed], 2), lo = r(gm[contrast == "registered" & group == "lo", G_mean_keyed_reversed], 2)),
                           congruence = list(groups = r(cg$mean_congruence_groups[2], 2), random = r(cg$mean_congruence_random_splits[2], 2))),
    background = list(r2 = .04, r2_background = .018, r2_effort = .019,
                      secs_or = r(exp(sr[model == "tells apart (log odds)" & term == "log_secs_z", estimate]), 2),
                      practice_pct = setNames(as.list(r(bg[variable == "practice", pct_tell_apart], 0)), bg[variable == "practice", level]),
                      r_best8 = r(scc[score == "best8" & with == "separation", r], 2), r_trial10 = r(scc[score == "trial10" & with == "separation", r], 2),
                      thirds = lapply(seq_len(nrow(sc)), function(i) list(third = as.character(sc$thirds[i]), best8 = r(sc$best8[i], 2))))),
  forms = lapply(seq_len(nrow(forms)), function(i) c(as.list(forms[i, .(form, n, n_reversed, omega_content = r(omega_content), share_agreement = r(share_agreement),
                                                                    alpha = r(alpha), alpha_told_apart = r(alpha_told_apart), alpha_not = r(alpha_not),
                                                                    mean = r(mean, 2), sd = r(sd, 2), mdc95 = r(mdc95, 2), d_practice = r(d_practice, 2),
                                                                    r_trial10 = r(r_trial10, 2), r_ius = r(r_ius, 2), r_self = r(r_self, 2), r_separation = r(r_separation, 2))]),
                                                     list(cfa = as.list(cfa[i, .(cfi = cfi.scaled, srmr, rmsea = rmsea.scaled, rmsea_hi = rmsea.ci.upper.scaled,
                                                                                 agreement_pct = agreement_pct_item_variance)])))),
  statements = lapply(seq_len(nrow(stm)), function(i) list(code = stm$source[i], text = stm$text[i], role = stm$role[i], written_for = stm$hypothesized[i],
                                                           reversed = stm$reversed[i], names_mode = stm$names_mode[i], loading = r(stm$G_keyed[i], 2),
                                                           cfa = r(if (stm$role[i] == "8") stm$cfa_8[i] else if (stm$role[i] == "10") stm$cfa_10[i] else stm$cfa_with_exploratory[i], 2),
                                                           agree = stm$agree_pct[i])),
  mode_free = list(final_u4 = as.list(mf[set == "final 8 (17)", .(omega = r(omega_content), alpha = r(alpha))]),
                   fully = as.list(mf[set == "final 8, both -> N6 + sfc_31", .(omega = r(omega_content), alpha = r(alpha), alpha_not = r(alpha_not))]),
                   forward_only = as.list(mf[set == "mode-free best 8, 0 reversed", .(omega = r(omega_content), agreement = r(share_agreement), alpha = r(alpha))]),
                   cfa_u4 = as.list(c20[model == "final 8 + RI", .(cfi = cfi.scaled, rmsea = rmsea.scaled)])),
  invariance = lapply(seq_len(nrow(inv)), function(i) as.list(inv[i, .(form, step, cfi = cfi.scaled, rmsea = rmsea.scaled, srmr)]))
)
writeLines(jsonlite::toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = NA), dest, useBytes = TRUE)
sf2a_log("wrote ", dest)
