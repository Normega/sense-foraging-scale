# 10: sensitivity analyses (Section 5, "Exploratory analysis"), each compared
# with the primary solution by Tucker congruence of the keyed content loadings:
#   1. the primary model in participants who tell the modes apart
#   2. the primary model without flagged participants
#   3. the primary model without the random-intercept factor (fitted in 05)
#   4. an EFA of the forward statements only (keyed, geomin, k factors)
# Also the step-crash group (as Study 1's translation-crash check, not registered).
# And, exploratory: the effect of presented position on statement means
# (Spearman r of answer with position) and on loadings (each statement's
# item-rest correlation with its factor when shown in the first vs second half).

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
A <- readRDS(sf2a_derived("sf2a"))
def <- A$definition; k <- SF2A$k_primary; X <- A$items_raw; M <- A$meta
P <- readRDS(sf2a_derived(sprintf("primary_k%02d", k)))
member <- P$member; fn <- P$fnames; L0 <- P$keyed
out <- list(); ns <- list()

refit <- function(keep, label) {
  sf2a_log(label, " (n = ", sum(keep), ")")
  ns[[label]] <<- sum(keep)
  # A subset can leave a statement with an empty category, which the ordinal
  # estimator cannot fit; that is recorded as the result, not a crash.
  m <- tryCatch(sf2a_fit_esem(X[keep, ], member, ri = TRUE), error = function(e) e)
  if (sf2a_failed(m)) {
    sf2a_log("  not estimable: ", sf2a_msg(m))
    return(data.table(factor = fn, matched = NA_character_, congruence = NA_real_))
  }
  sf2a_congruence(L0, sf2a_orient(sf2a_loadings(m, fn, def$items), member, def$reversed)$keyed)
}
out$tell_apart <- refit(M$tells_apart %in% TRUE, "1. participants who tell the modes apart")
out$without_flagged <- refit(!M$flagged, "2. without flagged participants")
out$without_step_crash <- refit(!M$step_crash, "   without the step-crash group")
if (!sf2a_failed(P$without_ri)) {
  out$without_random_intercept <- sf2a_congruence(L0, sf2a_orient(sf2a_loadings(P$without_ri, fn, def$items), member, def$reversed)$keyed)
  ns$without_random_intercept <- nrow(X)
}
fwd <- setdiff(def$items, def$reversed)
sf2a_log("4. EFA of the forward statements only")
m4 <- tryCatch(efa(data = as.data.frame(A$items_keyed[, fwd]), nfactors = k, ordered = TRUE, estimator = "WLSMV",
                   rotation = "geomin", missing = "pairwise", output = "lavaan"), error = function(e) e)
if (!inherits(m4, "error")) { out$forward_only <- sf2a_congruence(L0[fwd, ], inspect(m4, "std")$lambda); ns$forward_only <- nrow(X) }

res <- rbindlist(out, idcol = "analysis")
res[, n := unlist(ns)[match(analysis, names(ns))]]
fwrite(res, sf2a_out(sprintf("10_sensitivity_k%02d.csv", k)))
print(dcast(res, factor ~ analysis, value.var = "congruence"))

# Position effects. Positions were randomized, so a reliable trend is an order effect.
K <- A$items_keyed; POS <- A$positions[, def$items]
prim <- colnames(L0)[apply(abs(L0), 1, which.max)]
pe <- rbindlist(lapply(def$items, function(i) {
  others <- setdiff(def$items[prim == prim[def$items == i]], i)
  rest <- if (length(others)) rowMeans(K[, others, drop = FALSE], na.rm = TRUE) else rep(NA_real_, nrow(K))
  early <- POS[, i] <= 31
  ir <- function(s) suppressWarnings(cor(K[s, i], rest[s], use = "complete.obs"))
  data.table(item = i, r_position = round(suppressWarnings(cor(K[, i], POS[, i], use = "complete.obs", method = "spearman")), 3),
             item_rest_first_half = round(ir(early), 3), item_rest_second_half = round(ir(!early), 3))
}))
pe[, item_rest_diff := item_rest_second_half - item_rest_first_half]
fwrite(pe, sf2a_out("10_position_effects.csv"))
sf2a_log("position effects: median |r| = ", round(median(abs(pe$r_position), na.rm = TRUE), 3),
         "; median item-rest change, second half minus first = ", round(median(pe$item_rest_diff, na.rm = TRUE), 3))
