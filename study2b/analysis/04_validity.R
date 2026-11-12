# 04: H2-H3, convergent and discriminant validity (Section 5, model 2), and the
# predicted validity matrix.
#
# One model per validity measure: the H1 model plus that measure as latent factors of
# its items (subscales of one measure as correlated factors), correlated with the
# content factor S and uncorrelated with the agreement factor; 95% Wald CIs.
#   H2   the four lower bounds (FFMQ Observing, SBI, BFI-2 Aesthetic Sensitivity,
#        Brief MAIA-2 Noticing) exceed .30
#   H3a  their upper bounds are below .80 (Ronkko & Cho, 2022: "no problem")
#   H3b  five contrasts, each in a model holding both measures: the difference of the two
#        latent correlations (a defined parameter) has a 95% delta-method CI above 0;
#        H3b holds when at least four of the five do
#   H3c  both BIDR CIs lie within +/-.30
#   H3d  the multiple correlation of S with the four H2 neighbours (one model holding
#        all four), R = sqrt(r' Sigma^-1 r), 2,000-resample bootstrap 95% CI; holds when
#        its upper bound is below .80
#   H3e  BSSS-8's upper bound below .40; NCS-6, SPSQ-SF Sensory Discomfort and ASI-3
#        Physical Concerns within +/-.30
# Reported: every measure's latent correlation with the 8 (and the 10), the agreement
# factor's correlation with every measure (the same models with that correlation
# freed), the observed-score matrix, and the confirmatory cells without flagged
# participants. PREREG: single-item measures (life satisfaction, the past-week items,
# the engaged pause) enter as single-indicator factors with no residual variance, so the
# correlation is with the item's latent response. PREREG: contrasts 1, 2 and 4 come
# from the measure's own model, 3 and 5 from a model of just their two measures; H3d's
# CI is the bootstrap percentile interval, failed resamples left out and counted.

source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")
A <- readRDS(sf2b_derived("sf2b")); D <- A$defs; M <- A$meta; S <- A$scores
score_items <- sf2b_score_items(D)

DAT <- as.data.frame(A$scale_raw)
for (k in names(A$val_items)) if (k != "past_week") DAT <- cbind(DAT, A$val_items[[k]])
SINGLE <- c("lifesat", "crit_pause_days", "crit_notice_first", "crit_outdoors", "engaged_pause")
for (v in SINGLE) DAT[[v]] <- S[[v]]
rownames(DAT) <- NULL

# Contrast definitions added to the model that holds both measures.
contrast_lines <- function(scores) {
  h <- SF2B_H3B[first %in% scores & second %in% scores]
  if (!nrow(h)) return(character())
  sprintf("dif_%d := r_%s - r_%s", h$contrast, h$first, h$second)
}
TASKS <- list()
add_task <- function(...) TASKS[[length(TASKS) + 1L]] <<- list(...)
GROUPS <- c(lapply(SF2B_MEASURES, identity),
            list(pair3 = c("sbi_moment", "bfi_openmindedness"), pair5 = c("ffmq_observing", "samd_md")))
for (form in c("sf8", "sf10")) for (sample in c("all", "unflagged")) {
  if (form == "sf10" && sample == "unflagged") next
  for (g in names(GROUPS)) for (ri_free in c(FALSE, TRUE)) {
    if (sample == "unflagged" && ri_free) next
    if (g %in% c("pair3", "pair5") && ri_free) next
    add_task(form = form, sample = sample, group = g, scores = GROUPS[[g]], ri_free = ri_free, single = NULL)
  }
  for (v in SINGLE) for (ri_free in c(FALSE, TRUE)) {
    if (sample == "unflagged") next
    add_task(form = form, sample = sample, group = v, scores = character(), ri_free = ri_free, single = v)
  }
}
# H3d's model: the four H2 neighbours together.
add_task(form = "sf8", sample = "all", group = "h3d", scores = SF2B_H2, ri_free = FALSE, single = NULL)
add_task(form = "sf8", sample = "unflagged", group = "h3d", scores = SF2B_H2, ri_free = FALSE, single = NULL)

model_of <- function(t) {
  items <- D[[t$form]]
  syn <- c(sf2b_scale_syntax(items), sf2b_validity_syntax(t$scores, D, ri_free = t$ri_free), contrast_lines(t$scores))
  ord <- c(items, unlist(score_items[t$scores]))
  if (!is.null(t$single)) {
    v <- t$single
    syn <- c(syn, paste0("L_", v, " =~ 1*", v), paste0(v, " ~~ 0*", v), paste0("L_", v, " ~~ 1*L_", v),
             paste0("S ~~ r_", v, "*L_", v), paste0("RI ~~ ", if (t$ri_free) paste0("ri_", v, "*") else "0*", "L_", v))
    ord <- c(ord, v)
  }
  list(syn = syn, ord = unique(ord))
}
run_task <- function(t) {
  m <- model_of(t)
  d <- if (t$sample == "unflagged") DAT[!flagged, , drop = FALSE] else DAT
  f <- sf2b_try(sf2b_sem(m$syn, d, ordered = m$ord))
  if (sf2b_failed(f)) return(list(task = t, ok = FALSE, error = sf2b_msg(f)))
  bad <- sf2b_orientation(f, D[[t$form]], t$scores, D)
  if (length(bad)) return(list(task = t, ok = FALSE, error = paste("factor not oriented as keyed:", paste(bad, collapse = ", "))))
  pe <- parameterEstimates(f)
  ss <- standardizedSolution(f)
  r <- pe[grepl("^(r|dif)_", pe$label), c("label", "est", "se", "ci.lower", "ci.upper", "pvalue")]
  ri <- if (t$ri_free) { x <- ss[grepl("^ri_", ss$label), c("label", "est.std", "se", "ci.lower", "ci.upper", "pvalue")]
    names(x)[2] <- "est"; x } else NULL
  lv <- if (t$group == "h3d") lavInspect(f, "cor.lv") else NULL
  list(task = t, ok = TRUE, r = r, ri = ri, cor_lv = lv, fit = sf2b_fit_stats(f))
}
sf2b_log("fitting ", length(TASKS), " validity models on ", SF2B$cores, " workers")
RES <- sf2b_par_lapply(TASKS, run_task, data = list(DAT = DAT, D = D, score_items = score_items, flagged = M$flagged,
                                                    model_of = model_of, contrast_lines = contrast_lines))
failed <- Filter(function(r) !r$ok, RES)
if (length(failed)) for (r in failed) sf2b_log("FAILED: ", r$task$form, " ", r$task$sample, " ", r$task$group, if (r$task$ri_free) " (RI free)", ": ", r$error)

collect <- function(what) rbindlist(lapply(Filter(function(r) r$ok && !is.null(r[[what]]) && nrow(r[[what]]), RES), function(r)
  cbind(data.table(form = r$task$form, sample = r$task$sample, model = r$task$group, ri_free = r$task$ri_free), as.data.table(r[[what]]))))
R <- collect("r"); RIc <- collect("ri")
R[, score := sub("^r_", "", label)]
# The matrix: each score's correlation from its measure's own model (the primary
# estimate); the pair models supply only the two contrasts they exist for.
lat <- R[ri_free == FALSE & grepl("^r_", label) & !(model %in% c("pair3", "pair5", "h3d"))]
lat <- merge(lat, SF2B_PREDICTED, by = "score", all.x = TRUE, sort = FALSE)
lat[, (c("est", "se", "ci.lower", "ci.upper")) := lapply(.SD, r3), .SDcols = c("est", "se", "ci.lower", "ci.upper")]
fwrite(lat[order(form, sample, match(score, SF2B_PREDICTED$score))], sf2b_out("04_latent_validity.csv"))
RIc[, score := sub("^ri_", "", label)]
fwrite(RIc[, .(form, model, score, r_with_agreement = r3(est), lower = r3(ci.lower), upper = r3(ci.upper), p = signif(pvalue, 3))],
       sf2b_out("04_agreement_factor_validity.csv"))
fits <- rbindlist(lapply(Filter(function(r) r$ok, RES), function(r) c(list(form = r$task$form, sample = r$task$sample, model = r$task$group,
                                                                        ri_free = r$task$ri_free), r$fit)), fill = TRUE)
fwrite(fits, sf2b_out("04_model_fit.csv"))

# ---- H3d: the multiple correlation and its bootstrap -------------------------------------------
mult_r <- function(C) { r <- C["S", SF2B_H2]; Sig <- C[SF2B_H2, SF2B_H2]; sqrt(drop(t(r) %*% solve(Sig, r))) }
h3d_point <- function(sample) { x <- Filter(function(r) r$ok && r$task$group == "h3d" && r$task$sample == sample, RES)
  if (length(x)) mult_r(x[[1]]$cor_lv) else NA_real_ }
R_h3d <- h3d_point("all"); R_h3d_unfl <- h3d_point("unflagged")
B <- SF2B$boot_h3d
idx <- sf2b_boot_indices(nrow(DAT), B, SF2B$seed_boot)
t3d <- TASKS[[which(vapply(TASKS, function(t) t$group == "h3d" && t$sample == "all", TRUE))]]
m3d <- model_of(t3d)
sf2b_log("H3d bootstrap: ", B, " resamples")
boot <- unlist(sf2b_par_lapply(idx, function(ix) {
  f <- sf2b_try(sf2b_sem(m3d$syn, DAT[ix, , drop = FALSE], ordered = m3d$ord, fast = TRUE))
  if (sf2b_failed(f)) return(NA_real_)
  C <- lavInspect(f, "cor.lv"); r <- C["S", h2]; Sig <- C[h2, h2]
  v <- tryCatch(drop(t(r) %*% solve(Sig, r)), error = function(e) NA_real_)
  if (is.na(v) || v < 0) NA_real_ else sqrt(v)
}, data = list(m3d = m3d, DAT = DAT, h2 = SF2B_H2)))
fwrite(data.table(resample = seq_along(boot), R = boot), sf2b_out("04_h3d_bootstrap.csv"))
h3d_ci <- unname(quantile(boot, c(.025, .975), na.rm = TRUE))

# ---- decisions -------------------------------------------------------------------------------------
# Reported, not tested: Ronkko & Cho's (2022, Table 12) class of the 95% upper limit (the
# absolute lower limit when the estimate is negative).
rc_class <- function(est, lower, upper) { ul <- if (isTRUE(est < 0)) abs(lower) else upper
  if (is.null(ul) || is.na(ul)) NA_character_ else if (ul < .80) "no problem" else if (ul < .90) "marginal" else if (ul < 1) "moderate" else "severe" }
cell <- function(sc, fm = "sf8", sm = "all") lat[score == sc & form == fm & sample == sm]
dif_of <- function(k, fm = "sf8", sm = "all") R[label == paste0("dif_", k) & form == fm & sample == sm & ri_free == FALSE]
decide <- function(fm, sm) {
  h2 <- rbindlist(lapply(SF2B_H2, function(s) { x <- cell(s, fm, sm); data.table(score = s, r = x$est, lower = x$ci.lower, upper = x$ci.upper,
                                                                                   h2 = isTRUE(x$ci.lower > SF2B$h2_lower), h3a = isTRUE(x$ci.upper < SF2B$h3a_upper),
                                                                                   rc_class = rc_class(x$est, x$ci.lower, x$ci.upper)) }))
  h3b <- rbindlist(lapply(seq_len(nrow(SF2B_H3B)), function(k) { x <- dif_of(k, fm, sm); h <- SF2B_H3B[k]
    data.table(contrast = k, first = h$first, second = h$second, model = h$model,
               r_first = cell(h$first, fm, sm)$est, r_second = cell(h$second, fm, sm)$est,
               difference = r3(x$est), lower = r3(x$ci.lower), upper = r3(x$ci.upper), holds = isTRUE(x$ci.lower > 0)) }))
  h3c <- rbindlist(lapply(SF2B_H3C, function(s) { x <- cell(s, fm, sm)
    data.table(score = s, r = x$est, lower = x$ci.lower, upper = x$ci.upper, holds = isTRUE(x$ci.lower > -SF2B$h3c_bound && x$ci.upper < SF2B$h3c_bound)) }))
  h3e <- rbindlist(lapply(seq_len(nrow(SF2B_H3E)), function(k) { b <- SF2B_H3E[k]; x <- cell(b$score, fm, sm)
    data.table(score = b$score, r = x$est, lower = x$ci.lower, upper = x$ci.upper, bound_lower = b$lower, bound_upper = b$upper,
               holds = isTRUE(x$ci.lower > b$lower && x$ci.upper < b$upper)) }))
  list(h2 = h2, h3b = h3b, h3c = h3c, h3e = h3e,
       H2 = all(h2$h2), H3a = all(h2$h3a), H3b = sum(h3b$holds) >= SF2B$h3b_needed, H3b_contrasts_holding = sum(h3b$holds),
       H3c = all(h3c$holds), H3e = all(h3e$holds))
}
dec <- list(sf8 = decide("sf8", "all"), sf8_unflagged = decide("sf8", "unflagged"), sf10 = decide("sf10", "all"))
dec$sf8$h3d <- list(R = r3(R_h3d), lower = r3(h3d_ci[1]), upper = r3(h3d_ci[2]), resamples = B, resamples_ok = sum(!is.na(boot)),
                    holds = isTRUE(h3d_ci[2] < SF2B$h3d_upper), rc_class = rc_class(R_h3d, h3d_ci[1], h3d_ci[2]))
dec$sf8$H3d <- dec$sf8$h3d$holds
dec$sf8_unflagged$h3d <- list(R = r3(R_h3d_unfl), note = "point estimate only; the bootstrap is run on the full sample")
write_json(dec, sf2b_out("04_validity_decisions.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)

# ---- observed-score matrix (secondary) -------------------------------------------------------------------
obs_scores <- c(SF2B_PREDICTED$score, "bidr_sde_dichotomous", "bidr_im_dichotomous", "trial_av10")
obs <- rbindlist(lapply(c("sf8", "sf10", "trial_av10"), function(x) rbindlist(lapply(setdiff(obs_scores, x), function(s) {
  if (is.null(S[[s]])) return(NULL)
  r <- sf2b_r_ci(S[[x]], S[[s]]); data.table(score_of = x, with = s, n = r$n, r = r3(r$r), lower = r3(r$lower), upper = r3(r$upper), p = signif(r$p, 3)) }))))
fwrite(obs, sf2b_out("04_observed_validity.csv"))
saveRDS(list(results = RES, latent = lat, decisions = dec, boot = boot), sf2b_derived("validity"))
print(dec$sf8[c("H2", "H3a", "H3b", "H3b_contrasts_holding", "H3c", "H3d", "H3e")])
print(dec$sf8$h3b)
