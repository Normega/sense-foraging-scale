# 10: step 8, the revised form, if step 7 produced one (exploratory). Fitted on the
# HOLDOUT half only, by the H1 criteria and the per-factor rule; its omega and validity
# correlations are reported beside the fixed form's, in the same half. An internal
# check only: the fixed 20-item form stays the instrument of record.
# PREREG: "its validity correlations" are the observed-score Pearson correlations of
# each subscale with every validity score (the secondary analysis' matrix), in the
# holdout half's validity sample, for both forms.

source("study2/analysis/config.R")
source("study2/analysis/sf2_models.R")
A <- readRDS(sf2_derived("sf2")); S7 <- readRDS(sf2_derived("step7"))
D <- A$defs; M <- A$meta; F4 <- D$factors4
if (!nrow(S7$replacements)) {
  write_json(list(revised_form = FALSE, note = "step 7 selected no replacement: no revised form to check"),
             sf2_out("10_revised_form.json"), auto_unbox = TRUE, pretty = TRUE)
  sf2_log("no revised form")
} else {
hold <- S7$half[M$participant] == "holdout"
RAW <- cbind(A$sf_raw, A$cand_raw)
forms <- list(fixed = lapply(setNames(F4, F4), function(f) D$sf_items[D$sf_factor == f]), revised = S7$revised)
rev_all <- c(D$sf_rev, D$cand[reversed == TRUE, id])

fit_form <- function(form) {
  its <- unlist(form, use.names = FALSE); fo <- setNames(rep(names(form), lengths(form)), its)
  mem <- sf2_member(its, fo, F4)
  f <- sf2_try(sf2_sem(sf2_esem_syntax(its, F4, TRUE), as.data.frame(RAW[hold, its]), its, target = sf2_target(mem)))
  if (sf2_failed(f)) return(list(error = sf2_msg(f)))
  K <- sf2_orient(sf2_loadings(f, F4, its), mem, rev_all)$keyed
  rule <- sf2_factor_rule(K, mem); fs <- sf2_fit_stats(f); ri <- sf2_ri(f)
  list(h1 = sf2_h1_decision(fs, ri), fit = fs, ri_variance = ri$variance, rule = rule,
       loadings = data.table(item = its, factor = fo[its], r3(K)))
}
res <- lapply(forms, fit_form)

VS <- hold & M$reached_end_validity
vars <- SF2_VALIDITY_SCORES$score
VC <- rbindlist(lapply(names(forms), function(fm) rbindlist(lapply(F4, function(f) {
  sc <- sf2_score(RAW, forms[[fm]][[f]], rev_all)
  rbindlist(lapply(vars, function(v) { x <- sf2_r_ci(sc[VS], A$scores[[v]][VS])
    data.table(form = fm, factor = f, score = v, n = x$n, r = x$r, lower = x$lower, upper = x$upper) })) }))))
VCw <- dcast(VC, factor + score ~ form, value.var = "r")
VCw[, difference := revised - fixed]
fwrite(VC[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)], sf2_out("10_validity_by_form.csv"))

cmp <- rbindlist(lapply(names(res), function(fm) if (!is.null(res[[fm]]$error)) data.table(form = fm, error = res[[fm]]$error) else
  cbind(form = fm, h1 = res[[fm]]$h1, res[[fm]]$rule[, .(factor, confirmed, omega, scored_for_individuals)])), fill = TRUE)
fwrite(cmp, sf2_out("10_forms_holdout.csv"))
out <- list(revised_form = TRUE, n_holdout = sum(hold), n_holdout_validity = sum(VS), items = forms$revised,
            replacements = S7$replacements, fit = lapply(res, `[[`, "fit"), h1 = lapply(res, `[[`, "h1"),
            factors = cmp, validity_largest_differences = VCw[order(-abs(difference))][1:10])
write_json(out, sf2_out("10_revised_form.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null", dataframe = "rows")
print(cmp)
}
