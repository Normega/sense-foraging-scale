# 06: measurement invariance (Section 5, analysis 4; exploratory, specified in advance),
# on everyone analysable.
#
# Groups: posting; gender (men vs women); country group (per the Section 4 rule, as
# applied in 01; not analysed if it fails). Models: the four-factor CFA with the random
# intercept (20 items) and the two-factor Access + Valuing CFA with the random
# intercept (its 10 items). Steps (Wu & Estabrook 2016): configural; thresholds;
# + loadings; + intercepts (sf2_inv_syntax). Non-invariance: delta CFI < -.010 with
# delta RMSEA > .015 (scaled); then partial invariance item by item (sf2_inv_sequence).
# Sparse categories: a category with < 5 answers in any group is merged toward the
# centre for that item in all groups.
# ESEM invariance is reported when a CFA is not configurally acceptable (and always in
# simulated mode, so that the dry run exercises it). PREREG:
# "acceptable" = the H1 fit criteria (scaled CFI >= .95, RMSEA <= .06 with upper CI
# <= .08, SRMR <= .08). It is fitted as ESEM-within-CFA: the full-sample target-rotated
# ESEM's solution, with each factor's referent item (its highest own loading net of
# its largest cross-loading) keeping its cross-loadings fixed at the ESEM estimates in
# every group; everything else free. Its configural fit equals a multi-group ESEM's.
# Latent means (exploratory): posting and country, where the intercepts step holds.
# Over time: longitudinal invariance of Access + Valuing between baseline and follow-
# up (the same steps, each item's residuals correlated across occasions), follow-up
# sample.

source("study2/analysis/config.R")
source("study2/analysis/sf2_models.R")
A <- readRDS(sf2_derived("sf2"))
D <- A$defs; M <- A$meta; F4 <- D$factors4
sf20 <- D$sf_items; av10 <- sf20[D$sf_factor[sf20] %in% c("Access", "Valuing")]
fit_ok <- function(fs) isTRUE(fs$cfi >= .95 && fs$rmsea <= .06 && fs$rmsea_upper <= .08 && fs$srmr <= .08)

pattern_cfa <- function(its, fn) sf2_pattern_cfa(its, D$sf_factor, fn)
pattern_ewc <- function(its, fn) {
  m <- sf2_member(its, D$sf_factor, fn)
  e <- sf2_sem(sf2_esem_syntax(its, fn, TRUE), as.data.frame(A$sf_raw[, its]), its, target = sf2_target(m))
  Lu <- lavInspect(e, "est")$lambda[its, fn, drop = FALSE]; Ls <- sf2_loadings(e, fn, its)
  P <- matrix("free", length(its), length(fn), dimnames = list(its, fn))
  for (f in fn) {
    own <- its[m[, f]]
    score <- abs(Ls[own, f]) - apply(abs(Ls[own, setdiff(fn, f), drop = FALSE]), 1, max)
    ref <- own[which.max(score)]
    for (g in setdiff(fn, f)) P[ref, g] <- sprintf("%.6f", Lu[ref, g])
  }
  P
}

run_mg <- function(label, its, P, group, start0 = NULL) sf2_run_mg(A$sf_raw, its, P, group, label, start0)

groupings <- list(posting = M$posting, gender = M$gender2,
                  country = if (isTRUE(A$country_rule$analysed)) M$country_group else NULL)
groupings <- Filter(Negate(is.null), groupings)
models <- list(cfa20 = list(items = sf20, fn = F4), av10 = list(items = av10, fn = c("Access", "Valuing")))
runs <- list(); ewc_patterns <- list()
for (g in names(groupings)) for (mname in names(models)) {
  mm <- models[[mname]]; lab <- paste(mname, g, sep = " / ")
  runs[[lab]] <- run_mg(lab, mm$items, pattern_cfa(mm$items, mm$fn), groupings[[g]])
  conf <- runs[[lab]]$table[step == "configural"]
  # In simulated mode the ESEM path runs regardless, so the dry run proves it works.
  if (SF2$mode == "simulated" || !fit_ok(list(cfi = conf$cfi, rmsea = conf$rmsea, rmsea_upper = conf$rmsea_upper, srmr = conf$srmr))) {
    if (is.null(ewc_patterns[[mname]])) ewc_patterns[[mname]] <- pattern_ewc(mm$items, mm$fn)
    lab2 <- paste0(sub("cfa20", "esem20", sub("av10", "esem_av10", mname)), " / ", g)
    sf2_log("  configural CFA: CFI ", conf$cfi, ", RMSEA ", conf$rmsea, if (SF2$mode == "simulated") " (simulated mode: ESEM path run regardless)" else
            " (not acceptable)", ": ESEM-within-CFA")
    runs[[lab2]] <- run_mg(lab2, mm$items, ewc_patterns[[mname]], groupings[[g]], start0 = runs[[lab]]$fits$configural)
  }
}

# Over time: Access + Valuing, baseline vs follow-up.
P2 <- A$retest$pairs
if (nrow(P2) >= 50) {
  T1 <- A$retest$t1[, av10]; T2 <- A$retest$t2[, av10]
  mg <- sf2_merge_sparse(rbind(T1, T2), rep(1:2, each = nrow(T1)))
  ncat <- apply(mg$X, 2, function(x) length(unique(x[!is.na(x)])))
  dl <- as.data.frame(cbind(mg$X[seq_len(nrow(T1)), ], mg$X[nrow(T1) + seq_len(nrow(T2)), ]))
  names(dl) <- c(paste0(av10, "_t1"), paste0(av10, "_t2"))
  mem <- sf2_member(av10, D$sf_factor, c("Access", "Valuing"))
  build <- sf2_closure(function(step, free) sf2_long_syntax(mem, ncat, step, free), mem = mem, ncat = ncat)
  fit <- sf2_closure(function(syn, ...) sf2_sem(syn, dl, names(dl), auto.fix.first = FALSE, ...), dl = dl)
  sf2_log("longitudinal invariance: Access + Valuing, ", nrow(dl), " pairs")
  runs[["av10 / over time"]] <- sf2_inv_sequence(build, fit, av10, "av10 / over time")
  runs[["av10 / over time"]]$merges <- if (nrow(mg$log)) cbind(analysis = "av10 / over time", mg$log) else NULL
}

tab <- rbindlist(lapply(runs, `[[`, "table"), fill = TRUE)
fwrite(tab, sf2_out("06_invariance.csv"))
fwrite(rbindlist(lapply(runs, `[[`, "merges"), fill = TRUE), sf2_out("06_category_merges.csv"))

# Latent means where the intercepts step holds (exploratory): posting and country.
lm_tab <- rbindlist(lapply(names(runs), function(lab) {
  r <- runs[[lab]]; f <- r$fits$intercepts
  if (!grepl("posting|country", lab) || !isTRUE(r$completed) || is.null(f) || sf2_failed(f)) return(NULL)
  fn <- colnames(r$pattern)
  lam <- lavInspect(f, "est")[[1]]$lambda[rownames(r$pattern), fn, drop = FALSE]
  o <- sf2_orient(lam, sf2_member(rownames(r$pattern), D$sf_factor, fn), D$sf_rev)$orient
  lms <- sf2_latent_means(f, fn, o)
  cbind(analysis = lab, groups = paste(names(r$groups), collapse = " vs "), lms)
}), fill = TRUE)
if (nrow(lm_tab)) fwrite(lm_tab, sf2_out("06_latent_means.csv"))

# The highest step that holds: each step's final model (full or partial) converged and
# not non-invariant, in order, stopping at the first that fails.
held <- function(tb) {
  h <- NA_character_
  for (st in SF2_INV_STEPS) {
    rr <- tb[sub(" \\(partial\\)", "", step) == st]
    if (!nrow(rr)) break
    last <- rr[.N]
    if (!isTRUE(last$converged) || isTRUE(last$noninvariant)) break
    h <- if (grepl("partial", last$step)) paste(st, "(partial)") else st
  }
  h
}
summ <- rbindlist(lapply(names(runs), function(lab) { r <- runs[[lab]]
  data.table(analysis = lab, groups = paste(names(r$groups), unlist(r$groups), sep = " = ", collapse = ", "),
             highest_step_held = held(r$table),
             completed = r$completed, released = paste(unlist(Map(function(k, v) if (length(v)) paste0(k, ": ", paste(v, collapse = " ")), names(r$free), r$free)), collapse = "; "))
}))
write_json(list(summary = summ, country_rule = A$country_rule), sf2_out("06_invariance_summary.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null", dataframe = "rows")
saveRDS(lapply(runs, function(r) r[c("table", "free", "completed", "groups", "pattern")]), sf2_derived("invariance"))
print(summ)
