# 05: the primary model (Section 5, step 4) and decisions K and F.
#
# ESEM on the 58 UNRECODED statements (WLSMV, ordered, pairwise): SF2A$k_primary
# correlated content factors, partially specified target rotation (for five
# factors toward the hypothesized membership, non-member loadings targeted to 0;
# for any other number toward the membership in 04's geomin EFA), plus the
# random-intercept agreement factor: loadings 1 on all 58, variance estimated,
# uncorrelated with the content factors. The same model without the random
# intercept is fitted alongside (sensitivity 3).
#
# Decision F (Forage): Forage stands if the chosen solution has a factor on which
# at least three Forage statements have primary loadings >= .40 and more than half
# of the statements with primary loadings >= .40 are Forage statements.
# PREREG: a statement's primary factor is the content factor with its largest
# absolute keyed loading (the agreement factor is not a content factor), and a
# primary loading counts toward ">= .40" only in the keyed direction (a statement
# loading against its keying does not define the factor); here and in 06 and 09.

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
A <- readRDS(sf2a_derived("sf2a"))
def <- A$definition; X <- A$items_raw; k <- SF2A$k_primary
member <- sf2a_membership(def, k); fn <- colnames(member)

sf2a_log("ESEM with random intercept, k = ", k)
m_ri <- sf2a_fit_esem(X, member, ri = TRUE)
sf2a_log("ESEM without random intercept, k = ", k)
m_no <- sf2a_fit_esem(X, member, ri = FALSE)
stopifnot(!sf2a_failed(m_ri))

O <- sf2a_orient(sf2a_loadings(m_ri, fn, def$items), member, def$reversed)
Kor <- O$keyed
lab <- sf2a_label(Kor, def)
prim <- colnames(Kor)[apply(abs(Kor), 1, which.max)]
pload <- Kor[cbind(def$items, prim)]
sec <- apply(abs(Kor), 1, function(r) sort(r, decreasing = TRUE)[2])

# Factor scores (oriented). PREREG: "Factor scores use available statements";
# lavPredict's EBM scores are used, and any case it does not score is NA.
fs_raw <- tryCatch(lavPredict(m_ri), error = function(e) e)
FS <- matrix(NA_real_, nrow(X), length(fn), dimnames = list(rownames(X), fn))
if (!inherits(fs_raw, "error")) {
  idx <- lavInspect(m_ri, "case.idx")
  FS[idx, ] <- sweep(fs_raw[, fn, drop = FALSE], 2, O$orient, `*`)
  sf2a_log("factor scores for ", length(idx), " of ", nrow(X), " participants")
} else sf2a_log("factor scores failed: ", conditionMessage(fs_raw))

saveRDS(list(with_ri = m_ri, without_ri = m_no, member = member, fnames = fn, orient = O$orient, keyed = Kor,
             labels = lab, scores = FS), sf2a_derived(sprintf("primary_k%02d", k)))

fit_tab <- rbind(data.table(model = "ESEM + random intercept", as.data.table(sf2a_fit_stats(m_ri))),
                 data.table(model = "ESEM, no random intercept", as.data.table(sf2a_fit_stats(m_no))))
fwrite(fit_tab, sf2a_out(sprintf("05_primary_fit_k%02d.csv", k)))
ri <- sf2a_ri(m_ri)
fwrite(ri, sf2a_out(sprintf("05_agreement_factor_k%02d.csv", k)))

Lt <- data.table(item = def$items, source = def$source[def$items], hypothesized = def$factor[def$items],
                 keying = ifelse(def$items %in% def$reversed, "reverse", "forward"), round(Kor, 3),
                 RI_std = round(lavInspect(m_ri, "std")$lambda[def$items, "RI"], 3),
                 primary_factor = prim, primary = round(pload, 3), secondary = round(sec, 3))
fwrite(Lt, sf2a_out(sprintf("05_primary_loadings_keyed_k%02d.csv", k)))
phi <- lavInspect(m_ri, "std")$psi[fn, fn] * outer(O$orient, O$orient)
fwrite(data.table(factor = fn, round(phi, 3)), sf2a_out(sprintf("05_factor_correlations_k%02d.csv", k)))

# Decision K, criterion 1, in the primary model; and what each factor is.
facs <- rbindlist(lapply(fn, function(f) {
  on <- def$items[prim == f & pload >= SF2A$load_min]
  comp <- table(factor(def$factor[on], levels = def$factors))
  data.table(factor = f, matches = lab[[f]], n_primary_ge_40 = length(on),
             composition = paste(sprintf("%s %d", names(comp), comp)[comp > 0], collapse = ", "),
             n_reversed_ge_40 = sum(on %in% def$reversed))
}))
fwrite(facs, sf2a_out(sprintf("05_factors_k%02d.csv", k)))

# Decision F.
forage <- def$items[def$factor[def$items] == "Forage"]
dec_f <- rbindlist(lapply(fn, function(f) {
  on <- def$items[prim == f & pload >= SF2A$load_min]
  n_forage <- sum(on %in% forage)
  data.table(factor = f, n_primary_ge_40 = length(on), n_forage_ge_40 = n_forage,
             forage_share = if (length(on)) round(n_forage / length(on), 3) else NA_real_,
             forage_stands_here = n_forage >= 3 && n_forage > length(on) / 2)
}))
fwrite(dec_f, sf2a_out(sprintf("05_decision_F_k%02d.csv", k)))
joins <- data.table(item = forage, source = def$source[forage], primary_factor = prim[match(forage, def$items)],
                    matches = lab[prim[match(forage, def$items)]], primary = round(pload[match(forage, def$items)], 3))
fwrite(joins, sf2a_out(sprintf("05_forage_statements_k%02d.csv", k)))

sf2a_log("criterion 1 (every factor >= 3 statements at >= .40): ", all(facs$n_primary_ge_40 >= 3),
         "; agreement variance ", ri$variance, " (p ", ri$p, "); Forage stands: ", any(dec_f$forage_stands_here))
print(fit_tab); print(facs); print(dec_f)
