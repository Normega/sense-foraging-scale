# 12: step 10, the 14-item Access + Valuing form (exploratory), full analysable sample.
# Access 8 (9, 7, 6, 2, 8, 14, 36R, 35R) and Valuing 6 (18, 12, 25, 27, 28, 44R); items
# 6, 2, 14, 25 and 28 are the re-run Study 1 items sfc_21-25 in the candidate block.
#   - the two-factor CFA with the random intercept, judged by the H1 fit criteria;
#   - ordinal omega per factor (this model's loadings, sf2_omega);
#   - known groups (H6's model and d) and observed validity correlations, beside the
#     fixed form's Access and Valuing subscales (and their 10-item combination);
#   - the five re-run items' loadings and means against Study 1's (config.R). PREREG:
#     Study 2's loadings are taken from this CFA and, beside it, from step 7's 39-item
#     ESEM (selection half), the model closest to Study 1's.

source("study2/analysis/config.R")
source("study2/analysis/sf2_models.R")
A <- readRDS(sf2_derived("sf2")); S7 <- readRDS(sf2_derived("step7"))
D <- A$defs; M <- A$meta; S <- A$scores
RAW <- cbind(A$sf_raw, A$cand_raw); rev_all <- c(D$sf_rev, D$cand[reversed == TRUE, id])
its <- unlist(SF2_AV14, use.names = FALSE); fo <- setNames(rep(names(SF2_AV14), lengths(SF2_AV14)), its)
mem <- sf2_member(its, fo, names(SF2_AV14))
f <- sf2_sem(sf2_cfa_syntax(mem, ri = TRUE), as.data.frame(RAW[, its]), its)
K <- sf2_orient(sf2_loadings(f, names(SF2_AV14), its), mem, rev_all)$keyed
fs <- sf2_fit_stats(f); ri <- sf2_ri(f); rule <- sf2_factor_rule(K, mem)
phi <- lavInspect(f, "std")$psi["Access", "Valuing"] * prod(sf2_orient(sf2_loadings(f, names(SF2_AV14), its), mem, rev_all)$orient)

sc14 <- data.table(Access14 = sf2_score(RAW, SF2_AV14$Access, rev_all), Valuing14 = sf2_score(RAW, SF2_AV14$Valuing, rev_all),
                   AccessValuing14 = sf2_score(RAW, its, rev_all))
cmp <- cbind(S[, .(Access, Valuing, AccessValuing)], sc14)
grp <- ifelse(M$practice_group %in% "dedicated", 1L, ifelse(M$practice_group %in% "never", 0L, NA_integer_))
KG <- rbindlist(lapply(names(cmp), function(s) cbind(score = s, sf2_kg(cmp[[s]], grp, data.table(posting = M$posting_contemplative)))))
VS <- M$reached_end_validity
VC <- rbindlist(lapply(names(cmp), function(s) rbindlist(lapply(SF2_VALIDITY_SCORES$score, function(v) {
  x <- sf2_r_ci(cmp[[s]][VS], S[[v]][VS]); data.table(score = s, validity = v, n = x$n, r = x$r, lower = x$lower, upper = x$upper) }))))
VCw <- dcast(VC, validity ~ score, value.var = "r")

# Re-run items: loadings and means against Study 1.
RR <- copy(SF2_STUDY1_RERUN)
RR[, `:=`(study2_loading_av14_cfa = K[cbind(candidate, factor)],
          study2_loading_step7_esem = S7$loadings[cbind(candidate, factor)],
          study2_mean = colMeans(RAW[, candidate], na.rm = TRUE), study2_n = colSums(!is.na(RAW[, candidate])))]
RR[, `:=`(loading_difference = study2_loading_step7_esem - study1_loading, mean_difference = study2_mean - study1_mean)]

out <- list(fit = fs, decision = sf2_h1_decision(fs, ri), ri = ri, latent_r_access_valuing = phi, factors = rule,
            loadings = data.table(item = its, factor = fo[its], r3(K)), known_groups = KG, rerun_items = RR)
fwrite(VC[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)], sf2_out("12_av14_validity.csv"))
fwrite(KG[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)], sf2_out("12_av14_known_groups.csv"))
fwrite(RR[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], sf2_out("12_av14_rerun_items.csv"))
write_json(out, sf2_out("12_av14.json"), auto_unbox = TRUE, pretty = TRUE, digits = 4, na = "null", dataframe = "rows")
print(unlist(fs[c("cfi.scaled", "rmsea.scaled", "srmr")])); print(rule); print(KG[, .(score, d = round(d, 3), p = signif(p, 3))]); print(RR)
