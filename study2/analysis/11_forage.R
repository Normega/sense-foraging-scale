# 11: step 9, Forage (exploratory), on the full analysable sample.
#   (a) the five Forage items' alpha (keyed, Pearson), mean inter-item correlation, and
#       the Forage score's correlation with the agreement index (and, because two of
#       the index's pairs are Forage items, with the index without them);
#   (b) whether they form a factor of their own: an ESEM with five content factors
#       (Forage targeted) plus the random intercept. PREREG: the text does not say which
#       items; the primary model is the 20 short-form items + the 5 Forage items, and
#       the model with all 44 items (step 7's 39 + Forage) is reported beside it. "Its
#       own factor" = Forage meets H1's per-factor rule (4 of 5 keyed loadings >= .40,
#       none loading more strongly elsewhere); each item's loadings on Valuing and
#       Access are reported, as Study 1's discovery items loaded with Valuing;
#   (c) the Forage score's correlations with the past-week pause item, practice
#       frequency (codes 0-4, never-practitioners 0), BFI-2 Intellectual Curiosity
#       (general curiosity, Forage's neighbour), the rest of BFI-2 Open-Mindedness
#       (PREREG: its other 8 items, Aesthetic Sensitivity + Creative Imagination, scored
#       as their mean) and SBI Savoring the Moment (Pearson, 95% CI; validity measures in
#       the validity sample).

source("study2/analysis/config.R")
source("study2/analysis/sf2_models.R")
A <- readRDS(sf2_derived("sf2"))
D <- A$defs; M <- A$meta; S <- A$scores
FOR <- D$cand[factor == "Forage", id]; rev_all <- c(D$sf_rev, D$cand[reversed == TRUE, id])
RAW <- cbind(A$sf_raw, A$cand_raw)
K <- RAW[, FOR]; K[, intersect(FOR, rev_all)] <- 7 - K[, intersect(FOR, rev_all)]
has <- rowSums(!is.na(A$cand_raw)) > 0
C <- cor(K[has, ], use = "pair"); p <- ncol(C)
fs <- sf2_score(RAW, FOR, rev_all)
a <- list(n = sum(has), alpha = (p / (p - 1)) * (1 - p / sum(C)), mean_inter_item_r = mean(C[upper.tri(C)]),
          r_agreement = sf2_r_ci(fs, S$agreement), r_agreement_without_forage_pairs = sf2_r_ci(fs, S$agreement_noforage),
          item_correlations = as.data.table(round(C, 3), keep.rownames = "item"))

F5 <- c(D$factors4, "Forage")
factor_of <- c(D$sf_factor, setNames(D$cand$factor, D$cand$id))
esem5 <- function(its) {
  mem <- sf2_member(its, factor_of, F5)
  f <- sf2_try(sf2_sem(sf2_esem_syntax(its, F5, TRUE), as.data.frame(RAW[, its]), its, target = sf2_target(mem)))
  if (sf2_failed(f)) return(list(error = sf2_msg(f)))
  O <- sf2_orient(sf2_loadings(f, F5, its), mem, rev_all); Kl <- O$keyed
  rule <- sf2_factor_rule(Kl, mem)
  phi <- lavInspect(f, "std")$psi[F5, F5] * outer(O$orient, O$orient)
  list(n_items = length(its), fit = sf2_fit_stats(f), forage_rule = rule[factor == "Forage"],
       forage_items = data.table(item = FOR, r3(Kl[FOR, , drop = FALSE]),
                                 strongest = F5[apply(abs(Kl[FOR, , drop = FALSE]), 1, which.max)]),
       forage_correlations = as.list(r3(phi["Forage", D$factors4])))
}
items25 <- c(D$sf_items, FOR)
items44 <- c(D$sf_items, D$cand$id)
b <- list(primary_25_items = esem5(items25), sensitivity_44_items = esem5(items44))
b$forms_own_factor <- isTRUE(b$primary_25_items$forage_rule$confirmed)

VS <- M$reached_end_validity
bv <- D$validity$bfi2
om_rest <- sf2_score(A$val_items$bfi2, c(bv$scores$bfi_aesthetic, bv$scores$bfi_creative), bv$reverse, bv$range, "mean")
cc <- list(crit_pause_days = sf2_r_ci(fs[VS], S$crit_pause_days[VS]), practice_frequency = sf2_r_ci(fs, M$freq_code),
           bfi_intellectual_curiosity = sf2_r_ci(fs[VS], S$bfi_intellectual[VS]),
           bfi_openmindedness_rest = sf2_r_ci(fs[VS], om_rest[VS]), sbi_moment = sf2_r_ci(fs[VS], S$sbi_moment[VS]))
out <- list(items = FOR, a = a, b = b, c = cc)
write_json(out, sf2_out("11_forage.json"), auto_unbox = TRUE, pretty = TRUE, digits = 4, na = "null", dataframe = "rows")
str(a[1:3]); print(b$primary_25_items$forage_items); print(b$primary_25_items$forage_rule)
sf2_log("Forage forms its own factor: ", b$forms_own_factor)
