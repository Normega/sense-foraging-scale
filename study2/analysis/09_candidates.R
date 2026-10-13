# 09: step 7, candidate items (Section 5; exploratory, selection rules fixed).
#
# The analysable sample is split at random into a selection half and a holdout half,
# stratified by posting, with a recorded seed (SF2$seed_split). On the selection half:
# a target-rotated ESEM of 39 items (the 20 short-form items, the 14 new non-Forage
# candidates, the 5 re-run Study 1 items), four content factors, each candidate
# targeted to its a-priori factor, the random intercept on all 39 unrecoded items.
# Sparse categories are merged within the half (fewer than 5 answers -> toward the
# centre). A candidate qualifies as a replacement on its intended factor when:
#   1. its keyed primary loading is >= .50 and no other loading is >= .30 (absolute);
#   2. swapping it for the factor's weakest current item (lowest keyed loading)
#      raises the factor's model-based omega (sf2_omega, from this model's loadings);
#   3. the factor keeps at least one item of each keying direction;
#   4. its "Prefer not to answer" rate is <= 5% (PREREG: in the selection half).
# At most two items per factor are replaced. PREREG: the text does not say which when
# more qualify, or how a second replacement is judged. Procedure: candidates are tried
# in order of their keyed primary loading; the first that meets 1-4 against the current
# set replaces its weakest item; then the same again against the revised set (its
# weakest item, its omega), once more at most.
# Bootstrap stability: SF2$boot_step7 resamples of the selection half (stratified by
# posting), each repeating merge, fit and selection; each replacement's selection rate
# is reported, and one selected in fewer than 60% of resamples is marked unstable.
# Letting be and wording (Section 5, H1's secondary test, "with the candidates"): this
# model plus a negative-wording factor on every reversed item among the 39; the same
# rule (Letting be's keyed loadings on 4+ of its 5 short-form items < .30 -> not
# separable), and the loadings of all its items, including the four new positive ones.
# PREREG: run on the selection half, as "step 7's model" is.

source("study2/analysis/config.R")
source("study2/analysis/sf2_models.R")
A <- readRDS(sf2_derived("sf2"))
D <- A$defs; M <- A$meta; F4 <- D$factors4
new14 <- D$cand[factor != "Forage" & is.na(study1_item), id]; rerun5 <- D$cand[!is.na(study1_item), id]
stopifnot(length(new14) == 14, length(rerun5) == 5)
items39 <- c(D$sf_items, new14, rerun5)
factor_of <- c(D$sf_factor, setNames(D$cand$factor, D$cand$id))
rev39 <- intersect(items39, c(D$sf_rev, D$cand[reversed == TRUE, id]))
member <- sf2_member(items39, factor_of, F4); tg <- sf2_target(member)
X <- cbind(A$sf_raw, A$cand_raw)[, items39]

set.seed(SF2$seed_split)
half <- rep("holdout", nrow(M))
for (p in sort(unique(M$posting))) { ix <- which(M$posting == p); half[ix[sample.int(length(ix), floor(length(ix) / 2))]] <- "selection" }
sel <- half == "selection"
sf2_log("split: selection ", sum(sel), ", holdout ", sum(!sel), " (stratified by posting, seed ", SF2$seed_split, ")")

pna_rate <- function(rows) {
  # "Prefer not to answer" rate per candidate among those in `rows` who submitted the candidates
  sub <- A$meta$participant[rows]; has <- rowSums(!is.na(A$cand_raw[sub, , drop = FALSE])) > 0
  colMeans(is.na(A$cand_raw[sub[has], , drop = FALSE]))
}

fit39 <- function(Xh, fast = FALSE, extra = NULL) {
  mg <- sf2_merge_sparse(Xh)
  f <- sf2_try(sf2_sem(c(sf2_esem_syntax(items39, F4, TRUE), extra), as.data.frame(mg$X), items39, target = tg, fast = fast))
  list(fit = f, merges = mg$log)
}
keyed_loadings <- function(f) sf2_orient(sf2_loadings(f, F4, items39), member, rev39)$keyed

select <- function(K, pna) {
  reps <- list(); cand_tab <- list()
  for (f in F4) {
    cur <- D$sf_items[D$sf_factor == f]
    cands <- items39[factor_of[items39] == f & !(items39 %in% D$sf_items)]
    prim <- K[cands, f]; other <- apply(abs(K[cands, setdiff(F4, f), drop = FALSE]), 1, max)
    r1 <- prim >= .50 & other < .30; r4 <- pna[cands] <= .05
    weak0 <- cur[which.min(K[cur, f])]
    eval0 <- vapply(cands, function(cc) { ns <- c(setdiff(cur, weak0), cc)
      c(sf2_omega(K[ns, f]) > sf2_omega(K[cur, f]), any(ns %in% rev39) && any(!ns %in% rev39)) }, logical(2))
    cand_tab[[f]] <- data.table(factor = f, candidate = cands, keyed_primary = prim, max_other = other, pna_rate = pna[cands],
                                rule1 = r1, rule4 = r4, weakest_current = weak0, rule2_vs_weakest = eval0[1, ], rule3_vs_weakest = eval0[2, ],
                                selected = FALSE, round = NA_integer_)
    used <- character()
    for (rnd in 1:2) {
      weak <- cur[which.min(K[cur, f])]; om0 <- sf2_omega(K[cur, f]); done <- FALSE
      for (cc in cands[order(-prim)]) {
        if (cc %in% used || !isTRUE(r1[cc] && r4[cc])) next
        ns <- c(setdiff(cur, weak), cc); om1 <- sf2_omega(K[ns, f])
        if (om1 > om0 && any(ns %in% rev39) && any(!ns %in% rev39)) {
          reps[[length(reps) + 1]] <- data.table(factor = f, round = rnd, candidate = cc, replaces = weak,
                                                 candidate_loading = K[cc, f], replaced_loading = K[weak, f], omega_before = om0, omega_after = om1)
          cand_tab[[f]][candidate == cc, `:=`(selected = TRUE, round = rnd)]
          cur <- ns; used <- c(used, cc); done <- TRUE; break
        }
      }
      if (!done) break
    }
  }
  list(replacements = rbindlist(reps), candidates = rbindlist(cand_tab))
}

# ---- the selection half ------------------------------------------------------------------------------
S7 <- fit39(X[sel, ])
if (sf2_failed(S7$fit)) stop("step 7 model did not fit: ", sf2_msg(S7$fit))
K <- keyed_loadings(S7$fit)
pna_sel <- pna_rate(sel)
res <- select(K, pna_sel)
fs <- sf2_fit_stats(S7$fit)
LT <- data.table(item = items39, factor = factor_of[items39], keying = ifelse(items39 %in% rev39, "reverse", "forward"),
                 set = ifelse(items39 %in% D$sf_items, "short form", ifelse(items39 %in% rerun5, "Study 1 re-run", "new")), r3(K),
                 RI_std = r3(lavInspect(S7$fit, "std")$lambda[items39, "RI"]))
setnames(LT, F4, paste0(F4, "_keyed"))
fwrite(LT, sf2_out("09_selection_loadings.csv"))
fwrite(res$candidates[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], sf2_out("09_candidates.csv"))
fwrite(res$replacements[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], sf2_out("09_replacements.csv"))
fwrite(data.table(t(unlist(fs))), sf2_out("09_selection_fit.csv"))
if (nrow(S7$merges)) fwrite(S7$merges, sf2_out("09_category_merges.csv"))
sf2_log("selection half: ", nrow(res$replacements), " replacement(s)")

# ---- bootstrap stability ----------------------------------------------------------------------------------
Xs <- X[sel, ]; post_s <- M$posting[sel]
B <- SF2$boot_step7
idx <- sf2_boot_indices(nrow(Xs), B, SF2$seed_boot_step7, strata = post_s)
pna_s <- function(ix) { sub <- M$participant[sel][ix]; has <- rowSums(!is.na(A$cand_raw[sub, , drop = FALSE])) > 0
  colMeans(is.na(A$cand_raw[sub[has], , drop = FALSE])) }
boot_fun <- function(ix) {
  r <- fit39(Xs[ix, ], fast = TRUE)
  if (sf2_failed(r$fit)) return(NULL)
  sl <- select(keyed_loadings(r$fit), pna_s(ix))$replacements
  if (nrow(sl)) paste(sl$factor, sl$candidate) else character()
}
t0 <- Sys.time()
boots <- sf2_par_lapply(idx, boot_fun, data = list(fit39 = sf2_closure(fit39, items39 = items39, F4 = F4, tg = tg),
  keyed_loadings = sf2_closure(keyed_loadings, F4 = F4, items39 = items39, member = member, rev39 = rev39),
  select = sf2_closure(select, F4 = F4, D = list(sf_items = D$sf_items, sf_factor = D$sf_factor), items39 = items39, factor_of = factor_of, rev39 = rev39),
  pna_s = sf2_closure(pna_s, M = M[, .(participant)], sel = sel, A = list(cand_raw = A$cand_raw)), Xs = Xs))
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
ok <- !vapply(boots, is.null, TRUE)
picks <- table(unlist(boots[ok]))
rates <- data.table(replacement = names(picks), times_selected = as.integer(picks))[, rate := times_selected / sum(ok)]
rates[, c("factor", "candidate") := tstrsplit(replacement, " ")]
rates[, selected_in_half := paste(factor, candidate) %in% paste(res$replacements$factor, res$replacements$candidate)]
rates[, unstable := selected_in_half & rate < .60]
setorder(rates, -selected_in_half, -rate)
fwrite(rates[, .(factor, candidate, selected_in_half, times_selected, rate = round(rate, 3), unstable)], sf2_out("09_bootstrap_rates.csv"))
reps <- res$replacements
if (nrow(reps)) { reps[, rate := rates$rate[match(paste(factor, candidate), paste(rates$factor, rates$candidate))]]
  reps[is.na(rate), rate := 0]; reps[, unstable := rate < .60] }

# ---- Letting be and wording, with the candidates ------------------------------------------------------------
wfit <- fit39(X[sel, ], extra = c(paste0("W =~ NA*", rev39[1], " + ", paste(rev39[-1], collapse = " + ")), "W ~~ 1*W",
                                  paste0("W ~~ ", paste0("0*", c(F4, "RI"), collapse = " + "))))
wording <- if (sf2_failed(wfit$fit)) list(error = sf2_msg(wfit$fit)) else {
  Kw <- keyed_loadings(wfit$fit); lb_sf <- D$sf_items[D$sf_factor == "LettingBe"]; lb_all <- items39[factor_of[items39] == "LettingBe"]
  list(fit = sf2_fit_stats(wfit$fit), letting_be_sf_items_below_30 = sum(Kw[lb_sf, "LettingBe"] < .30),
       letting_be_separable = sum(Kw[lb_sf, "LettingBe"] < .30) < 4,
       letting_be_loadings = as.list(r3(Kw[lb_all, "LettingBe"])), wording_loadings = as.list(r3(lavInspect(wfit$fit, "std")$lambda[rev39, "W"])))
}

revised <- lapply(setNames(F4, F4), function(f) { cur <- D$sf_items[D$sf_factor == f]
  r <- reps[factor == f]; for (k in seq_len(nrow(r))) cur <- c(setdiff(cur, r$replaces[k]), r$candidate[k]); cur })
out <- list(n_selection = sum(sel), n_holdout = sum(!sel), seed = SF2$seed_split, fit = fs,
            replacements = reps, revised_form = revised, any_replacement = nrow(reps) > 0,
            bootstrap = list(resamples = B, failed = sum(!ok), seconds = round(secs, 1),
                             core_seconds_per_resample = round(secs / B * min(SF2$cores, B), 2), cores = SF2$cores),
            wording_with_candidates = wording)
write_json(out, sf2_out("09_step7.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null", dataframe = "rows")
saveRDS(list(half = setNames(half, M$participant), revised = revised, replacements = reps, items39 = items39,
             loadings = K, factor_of = factor_of, rev39 = rev39), sf2_derived("step7"))
print(res$replacements); print(rates)
sf2_log("step 7 bootstrap: ", B, " resamples in ", round(secs), " s (", out$bootstrap$core_seconds_per_resample, " core-seconds each), ",
        sum(!ok), " failed")
