# 03: H1, structure (Section 5, analysis 1), on everyone analysable.
#
# Primary model: ESEM of the 20 UNRECODED short-form items, four correlated content
# factors, oblique partially specified target rotation (rotation = "pst", 30 random
# starts; each item targeted to its Section 4 factor, every other loading to 0), plus
# the random intercept. Decisions:
#   H1 supported: scaled CFI >= .95, RMSEA <= .06 (upper 90% CI <= .08), SRMR <= .08,
#     random-intercept variance > 0 (p < .05); partial: see sf2_h1_decision.
#   a factor is confirmed: >= 4 of its 5 items load >= .40 in their keyed direction and
#     none of its items loads more strongly on another factor;
#   scored for individuals: confirmed and model-based ordinal omega >= .70 (content
#     loadings only).
# Also reported: Model A (simple-structure CFA + RI) with its modification indices
# (nothing is re-specified on them), the ESEM without the RI, a one-factor model + RI.
# Pre-specified secondary: the primary ESEM plus a negative-wording factor (free
# loadings on the 8 reversed items, orthogonal to everything); Letting be is "not
# separable from reversed wording" if its keyed loadings on 4 or more of its 5 items
# fall below .30 in that model.

source("study2/analysis/config.R")
source("study2/analysis/sf2_models.R")
A <- readRDS(sf2_derived("sf2"))
D <- A$defs; items <- D$sf_items; F4 <- D$factors4
X <- as.data.frame(A$sf_raw)
member <- sf2_member(items, D$sf_factor, F4); tg <- sf2_target(member)

sf2_log("H1 primary: ESEM + random intercept, N = ", nrow(X))
fits <- list()
fits$primary <- sf2_sem(sf2_esem_syntax(items, F4, ri = TRUE), X, items, target = tg)
fits$esem_no_ri <- sf2_sem(sf2_esem_syntax(items, F4, ri = FALSE), X, items, target = tg)
fits$model_a <- sf2_sem(sf2_cfa_syntax(member, ri = TRUE), X, items)
one <- c(paste0("G =~ NA*", items[1], " + ", paste(items[-1], collapse = " + ")), "G ~~ 1*G", sf2_ri_syntax(items, "G"))
fits$one_factor <- sf2_sem(one, X, items)
rev_sf <- intersect(items, D$sf_rev)
wording <- c(sf2_esem_syntax(items, F4, ri = TRUE),
             paste0("W =~ NA*", rev_sf[1], " + ", paste(rev_sf[-1], collapse = " + ")), "W ~~ 1*W",
             paste0("W ~~ ", paste0("0*", c(F4, "RI"), collapse = " + ")))
fits$wording <- sf2_try(sf2_sem(wording, X, items, target = tg))

fit_tab <- rbindlist(lapply(names(fits), function(m) {
  f <- fits[[m]]; ri <- if (!sf2_failed(f)) sf2_ri(f) else NULL
  c(list(model = m, converged = !sf2_failed(f), error = sf2_msg(f)), sf2_fit_stats(f),
    list(ri_variance = ri$variance %||% NA, ri_se = ri$se %||% NA, ri_p = ri$p %||% NA,
         ri_pct_item_variance = ri$mean_pct_item_variance %||% NA))
}), fill = TRUE)
fwrite(fit_tab, sf2_out("03_fit.csv"))

# Primary model: loadings keyed and oriented, factor correlations, the per-factor rule.
P <- fits$primary
L <- sf2_loadings(P, F4, items); O <- sf2_orient(L, member, D$sf_rev)
ri_std <- lavInspect(P, "std")$lambda[items, "RI"]
load_tab <- data.table(item = items, factor = D$sf_factor[items], keying = ifelse(items %in% D$sf_rev, "reverse", "forward"),
                       r3(O$keyed), RI_std = r3(ri_std))
setnames(load_tab, F4, paste0(F4, "_keyed"))
fwrite(load_tab, sf2_out("03_h1_loadings.csv"))
phi <- lavInspect(P, "std")$psi[F4, F4] * outer(O$orient, O$orient)
fwrite(data.table(factor = F4, r3(phi)), sf2_out("03_h1_factor_correlations.csv"))
rule <- sf2_factor_rule(O$keyed, member)
fwrite(rule, sf2_out("03_h1_factor_rule.csv"))
ri <- sf2_ri(P); fs <- sf2_fit_stats(P)
decision <- list(h1 = sf2_h1_decision(fs, ri), fit = fs, ri = ri,
                 factors = rule[, .(factor, confirmed, omega, scored_for_individuals)])

# Modification indices for Model A: reported, nothing re-specified on them.
mi <- sf2_try(modindices(fits$model_a, sort. = TRUE, maximum.number = 40))
if (!inherits(mi, "error")) fwrite(as.data.table(mi)[, .(lhs, op, rhs, mi = round(mi, 2), epc = round(epc, 3), sepc.all = round(sepc.all, 3))],
                                   sf2_out("03_modelA_modification_indices.csv"))

# Letting be and wording.
W <- fits$wording
if (!sf2_failed(W)) {
  Lw <- sf2_loadings(W, F4, items); Ow <- sf2_orient(Lw, member, D$sf_rev)
  wl <- lavInspect(W, "std")$lambda[items, "W"]
  lb <- items[D$sf_factor[items] == "LettingBe"]
  wtab <- data.table(item = items, factor = D$sf_factor[items], keying = ifelse(items %in% D$sf_rev, "reverse", "forward"),
                     r3(Ow$keyed), wording_std = r3(wl))
  setnames(wtab, F4, paste0(F4, "_keyed"))
  n_below <- sum(Ow$keyed[lb, "LettingBe"] < .30)
  decision$wording <- list(model_fit = sf2_fit_stats(W), letting_be_items_below_30 = n_below,
                           letting_be_separable = n_below < 4, letting_be_loadings = as.list(r3(Ow$keyed[lb, "LettingBe"])))
  fwrite(wtab, sf2_out("03_wording_loadings.csv"))
} else decision$wording <- list(error = sf2_msg(W))

# Sensitivity (Section 5, "Data exclusion": flagged cases are kept, and flagged "for
# sensitivity analyses"): the primary model without anyone flagged. PREREG: the text
# names no analyses; H1, H2-H4 (04), H6 (07) and H7 (08) are each repeated without them.
unfl <- !A$meta$flagged
Pu <- sf2_try(sf2_sem(sf2_esem_syntax(items, F4, TRUE), X[unfl, ], items, target = tg))
decision$sensitivity_unflagged <- if (sf2_failed(Pu)) list(error = sf2_msg(Pu)) else {
  Ou <- sf2_orient(sf2_loadings(Pu, F4, items), member, D$sf_rev); fsu <- sf2_fit_stats(Pu); riu <- sf2_ri(Pu)
  list(n = sum(unfl), h1 = sf2_h1_decision(fsu, riu), fit = fsu, ri_variance = riu$variance,
       factors = sf2_factor_rule(Ou$keyed, member)[, .(factor, confirmed, omega, scored_for_individuals)],
       congruence_with_primary = as.list(round(diag(psych::factor.congruence(Ou$keyed, O$keyed)), 3)))
}

saveRDS(list(fits = fits, member = member, orient = O$orient, rule = rule, decision = decision), sf2_derived("h1"))
write_json(decision, sf2_out("03_h1_decision.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
print(fit_tab[, .(model, cfi.scaled, rmsea.scaled, rmsea.ci.upper.scaled, srmr, ri_variance)])
print(rule)
sf2_log("H1: ", decision$h1, "; Letting be separable from wording: ", decision$wording$letting_be_separable %||% NA)
