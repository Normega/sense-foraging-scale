# 04: H2-H4, validity (Section 5, analysis 2), in the validity sample (participants
# who reached the end of the validity block).
#
# Primary analysis: the H1 model, with each validity score added IN TURN as an
# observed variable (a single-indicator latent) correlated with the four content
# factors and uncorrelated with the random intercept; latent-observed correlations
# with 95% CIs (delta method). Validity scores are z-scored first (correlations do not
# change; the ordinal and continuous parts are then on comparable scales).
#   H2   each listed correlation in the predicted direction, 95% CI lower bound > .30
#   H3a  each H2 correlation's 95% CI upper bound < .80
#   H3b  Valuing's correlation with SBI exceeds its correlation with BFI-2 Open-
#        Mindedness: bootstrap difference in the latent model, 95% percentile CI
#        excluding 0. PREREG: both measures are added to the same model (a difference
#        between two correlations needs both), correlated with each other.
#   H3c  every 95% CI for the correlations with BIDR SDE and IM lies within +-.30
#   H3d  Access / Valuing with its neighbour set added (one model per factor): R =
#        sqrt(r' S^-1 r), r = the factor's correlations with the neighbours, S theirs
#        with each other; supported when the upper bound of the bootstrap 95% CI < .80
#   H4   the listed correlations positive, one-tailed alpha .05, Holm across the four
# Secondary: Pearson correlations of the observed short-form scores with every
# validity score (the full validity matrix). Exploratory: the random intercept's
# correlation with every validity score (the same model with that correlation free).

source("study2/analysis/config.R")
source("study2/analysis/sf2_models.R")
A <- readRDS(sf2_derived("sf2"))
D <- A$defs; items <- D$sf_items; F4 <- D$factors4
VS <- A$meta$reached_end_validity
vars <- SF2_VALIDITY_SCORES$score
Zs <- scale(as.matrix(A$scores[VS, ..vars]))
dat <- cbind(as.data.frame(A$sf_raw[VS, ]), as.data.frame(Zs))
member <- sf2_member(items, D$sf_factor, F4); tg <- sf2_target(member)
sf2_log("validity sample: N = ", nrow(dat))

# Fit the H1 model with observed variables vs; return oriented factor-variable
# correlations (and, if se, their CIs), RI-variable correlations and the variables'
# correlations with each other.
lat <- sf2_closure(function(vs, ri = "zero", fast = FALSE, d = dat) {
  f <- sf2_try(sf2_sem(c(sf2_esem_syntax(items, F4, TRUE), sf2_observed_syntax(vs, F4, ri)),
                       d[, c(items, vs)], items, target = tg, fast = fast))
  if (sf2_failed(f)) return(list(ok = FALSE, error = sf2_msg(f)))
  o <- sf2_orient(sf2_loadings(f, F4, items), member, rev_items)$orient
  psi <- lavInspect(f, "std")$psi; lv <- paste0("L_", vs)
  r <- psi[F4, lv, drop = FALSE] * o
  out <- list(ok = TRUE, r = r, S = psi[lv, lv, drop = FALSE], ri = if (ri == "free") psi["RI", lv] else NULL, orient = o)
  if (!fast) {
    ss <- as.data.table(standardizedSolution(f))
    ss <- ss[op == "~~" & ((lhs %in% c(F4, "RI") & rhs %in% lv) | (rhs %in% c(F4, "RI") & lhs %in% lv))]
    ss[, `:=`(factor = ifelse(lhs %in% lv, rhs, lhs), lvar = ifelse(lhs %in% lv, lhs, rhs))]
    ss[, s := ifelse(factor == "RI", 1, o[factor])]
    out$table <- ss[, .(score = sub("^L_", "", lvar), factor, r = est.std * s, se,
                        lower = ifelse(s > 0, ci.lower, -ci.upper), upper = ifelse(s > 0, ci.upper, -ci.lower),
                        z = z * s, p = pvalue)]
  }
  out
}, dat = dat, items = items, F4 = F4, tg = tg, member = member, rev_items = D$sf_rev)

# ---- primary: each validity score in turn (and the exploratory RI-free version) ------------------
t0 <- Sys.time()
prim <- sf2_par_lapply(vars, function(v) list(zero = lat(v, "zero"), free = lat(v, "free")), data = list(lat = lat))
names(prim) <- vars
sf2_log("latent validity: ", length(vars), " scores x 2 models in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")
LV <- rbindlist(lapply(vars, function(v) { p <- prim[[v]]$zero
  if (!p$ok) return(data.table(score = v, factor = F4, error = p$error)); p$table[factor %in% F4] }), fill = TRUE)
LV[, p_one_tailed := pnorm(z, lower.tail = FALSE)]
pred <- melt(SF2_VALIDITY_SCORES, id.vars = c("score", "label"), variable.name = "factor", value.name = "predicted")
LV <- merge(LV, pred, by = c("score", "factor"), all.x = TRUE, sort = FALSE)
band_ok <- function(r, code) mapply(function(r, code) { if (is.na(code) || is.na(r)) return(NA)
  b <- SF2_BANDS[[code]]; if (b[3] == 0) abs(r) < b[2] else sign(r) == b[3] && abs(r) >= b[1] && abs(r) <= b[2] }, r, code)
LV[, consistent_with_prediction := band_ok(r, predicted)]
RI <- rbindlist(lapply(vars, function(v) { p <- prim[[v]]$free
  if (!p$ok) return(data.table(score = v, error = p$error)); p$table[factor == "RI"] }), fill = TRUE)
fwrite(LV[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)], sf2_out("04_latent_validity.csv"))
fwrite(RI[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)], sf2_out("04_random_intercept_validity.csv"))

# ---- secondary: observed scores (Pearson, 95% CI) ---------------------------------------------------------
SFS <- c("Access", "Understanding", "Valuing", "LettingBe", "AccessValuing")
OBS <- rbindlist(lapply(vars, function(v) rbindlist(lapply(SFS, function(f) {
  x <- sf2_r_ci(A$scores[[f]][VS], A$scores[[v]][VS]); data.table(score = v, factor = f, n = x$n, r = x$r, lower = x$lower, upper = x$upper, p = x$p) }))))
fwrite(OBS[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)], sf2_out("04_observed_validity.csv"))

# ---- decisions: H2, H3a, H3c, H4 ---------------------------------------------------------------------------
cell <- function(f, v) LV[factor == f & score == v]
h2 <- rbindlist(lapply(seq_len(nrow(SF2_H2)), function(k) { x <- cell(SF2_H2$factor[k], SF2_H2$score[k])
  data.table(factor = x$factor, score = x$score, r = x$r, lower = x$lower, upper = x$upper,
             h2_supported = isTRUE(x$r > 0 && x$lower > .30), h3a_supported = isTRUE(x$upper < .80)) }))
h3c <- LV[score %in% c("bidr_sde", "bidr_im"), .(factor, score, r, lower, upper, within = lower > -.30 & upper < .30)]
h4 <- rbindlist(lapply(seq_len(nrow(SF2_H4)), function(k) { x <- cell(SF2_H4$factor[k], SF2_H4$score[k])
  data.table(factor = x$factor, score = x$score, r = x$r, z = x$z, p_one_tailed = x$p_one_tailed) }))
h4[, p_holm := p.adjust(p_one_tailed, "holm")][, supported := r > 0 & p_holm < .05]

# ---- H3b and H3d: the joint models, and their bootstraps -------------------------------------------------------
stat_h3b <- function(x) x$r["Valuing", "L_sbi_moment"] - x$r["Valuing", "L_bfi_openmindedness"]
stat_h3d <- function(x, f) { r <- x$r[f, ]; sqrt(drop(t(r) %*% solve(x$S, r))) }
full_b <- lat(c("sbi_moment", "bfi_openmindedness"))
full_d <- lapply(names(SF2_H3D), function(f) lat(SF2_H3D[[f]])); names(full_d) <- names(SF2_H3D)
B <- SF2$boot_h3
idx <- sf2_boot_indices(nrow(dat), B, SF2$seed_boot_h3)
boot_fun <- function(ix) {
  db <- dat[ix, ]
  b <- lat(c("sbi_moment", "bfi_openmindedness"), fast = TRUE, d = db)
  d <- lapply(names(nb), function(f) lat(nb[[f]], fast = TRUE, d = db)); names(d) <- names(nb)
  c(h3b = if (b$ok) stat_h3b(b) else NA, sapply(names(nb), function(f) if (d[[f]]$ok) stat_h3d(d[[f]], f) else NA))
}
t0 <- Sys.time()
boots <- sf2_par_lapply(idx, boot_fun, data = list(lat = lat, dat = dat, nb = SF2_H3D, stat_h3b = stat_h3b, stat_h3d = stat_h3d))
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
BT <- as.data.table(do.call(rbind, boots)); setnames(BT, c("h3b", paste0("h3d_", names(SF2_H3D))))
fwrite(BT[, lapply(.SD, round, 4)], sf2_out("04_bootstrap_draws.csv"))
ci <- function(v) as.list(quantile(v, c(.025, .975), na.rm = TRUE, names = FALSE))
h3b_ci <- ci(BT$h3b)
h3b <- list(r_sbi = full_b$r["Valuing", "L_sbi_moment"], r_open = full_b$r["Valuing", "L_bfi_openmindedness"],
            difference = stat_h3b(full_b), lower = h3b_ci[[1]], upper = h3b_ci[[2]],
            resamples = B, failed_resamples = sum(is.na(BT$h3b)), supported = isTRUE(stat_h3b(full_b) > 0 && h3b_ci[[1]] > 0))
h3d <- rbindlist(lapply(names(SF2_H3D), function(f) { cc <- ci(BT[[paste0("h3d_", f)]])
  data.table(factor = f, neighbours = paste(SF2_H3D[[f]], collapse = " "), R = stat_h3d(full_d[[f]], f), lower = cc[[1]], upper = cc[[2]],
             resamples = B, failed_resamples = sum(is.na(BT[[paste0("h3d_", f)]])), supported = isTRUE(cc[[2]] < .80),
             r_with_neighbours = paste(round(full_d[[f]]$r[f, ], 3), collapse = " ")) }))

# Sensitivity: the confirmatory cells (H2, H3c, H4) refitted without flagged participants.
unfl <- !A$meta$flagged[VS]
sens_vars <- unique(c(SF2_H2$score, "bidr_sde", "bidr_im", SF2_H4$score))
sens <- rbindlist(sf2_par_lapply(sens_vars, function(v) { x <- lat(v, d = du); if (!x$ok) data.table(score = v, error = x$error) else x$table[factor %in% F4] },
                                 data = list(lat = lat, du = dat[unfl, ], F4 = F4)), fill = TRUE)
sens[, p_one_tailed := pnorm(z, lower.tail = FALSE)]
fwrite(sens[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)], sf2_out("04_latent_validity_unflagged.csv"))
cellu <- function(f, v) sens[factor == f & score == v]
sens_dec <- list(n = sum(unfl),
  h2 = rbindlist(lapply(seq_len(nrow(SF2_H2)), function(k) { x <- cellu(SF2_H2$factor[k], SF2_H2$score[k])
    data.table(factor = x$factor, score = x$score, r = x$r, lower = x$lower, upper = x$upper, h2_supported = isTRUE(x$r > 0 && x$lower > .30),
               h3a_supported = isTRUE(x$upper < .80)) })),
  h3c_supported = all(sens[score %in% c("bidr_sde", "bidr_im"), lower > -.30 & upper < .30]),
  h4 = { h <- rbindlist(lapply(seq_len(nrow(SF2_H4)), function(k) { x <- cellu(SF2_H4$factor[k], SF2_H4$score[k])
    data.table(factor = x$factor, score = x$score, r = x$r, p_one_tailed = x$p_one_tailed) }))
    h[, p_holm := p.adjust(p_one_tailed, "holm")][, supported := r > 0 & p_holm < .05]; h })

dec <- list(n = nrow(dat), sensitivity_unflagged = sens_dec, h2 = h2, h2_all_supported = all(h2$h2_supported), h3a = h2[, .(factor, score, upper, h3a_supported)],
            h3a_all_supported = all(h2$h3a_supported), h3b = h3b, h3c = h3c, h3c_supported = all(h3c$within), h3d = h3d, h4 = h4,
            bootstrap_seconds = round(secs, 1), seconds_per_resample = round(secs / B * min(SF2$cores, B), 2), cores = SF2$cores)
write_json(dec, sf2_out("04_validity_decisions.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null", dataframe = "rows")
saveRDS(list(primary = prim, full_b = full_b, full_d = full_d, boot = BT, decisions = dec), sf2_derived("validity"))
print(h2); print(h3c); print(h4); print(h3d); str(h3b)
sf2_log("bootstrap: ", B, " resamples in ", round(secs), " s on ", SF2$cores, " cores (",
        round(secs / B * min(SF2$cores, B), 2), " core-seconds per resample, 3 fits each)")
