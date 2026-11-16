# 21: is there a subset of the 13 statements that separates better from FFMQ Observing
# and the other near neighbours while staying a sound scale? EXPLORATORY, NOT REGISTERED.
# POST HOC.
#
# Not part of the filed registration (osf.io/d4a5u) and not in its list of exploratory
# analyses. Specified on 16 November 2026, after the confirmatory real-mode run and
# analyses 11-19 had been seen (14: no single statement drives the overlap; U4, sfc_31
# and N6 are weak; the reversed statements carry a negative specific association with
# Observing that disappears when the agreement factor may correlate with Observing; 15:
# sf6_44 in place of U4 breaks the H1 rule). Any subset chosen here is a hypothesis for a
# new sample, not a result. No registered cell, test or output is changed, and the
# registered 8 remains the scale.
#
# Design.
#   Pool: the 13 statements (the 8, A9 and sfc_24, and the unscored reversed N6, sfc_31,
#   sf6_44). Candidates: every subset of 5 to 10 statements with at least one reversed
#   and at least one unreversed statement (the agreement factor, loadings 1 on the
#   unrecoded statements, is identified only against content running both ways): 6,971.
#   Split: the 720 analysable participants split at random into a search half and a test
#   half (360 each), stratified by practice group, 20 splits drawn in order from one
#   seed; split 1 is the primary split. The whole search is repeated in every split.
#   Search (search half only):
#     screen: every candidate's point estimates computed from the half's pairwise
#       polychoric matrix with lavaan's own DWLS weights (lavCor), by a direct DWLS fit of
#       the same theta-parameterized model (residual variances 1, S variance 1, the
#       agreement factor's loadings 1 and variance free): the H1 model (content omega,
#       largest residual), 04's FFMQ model (five facets; r with Observing) and 04's H3d
#       model (the four H2 neighbours; their four r and the multiple R). Thresholds are
#       saturated and the weights diagonal, so these are WLSMV's point estimates exactly
#       (checked against lavaan on the registered 8 and on every candidate checked).
#     constraints, in lavaan (the registered code): the H1 rule as registered (sf2b_h1),
#       content omega >= .75, and no improper solution (10's improper_of()) in the H1
#       model, the FFMQ model or the H3d model; every factor oriented as keyed.
#     selection: candidates with screened omega >= .75 are walked in order of the
#       objective, each checked in lavaan, until the five best that meet every
#       constraint are found (at most 400 checked per objective).
#   Objectives (each minimised; each under both agreement models):
#     1 r with FFMQ Observing (04's FFMQ model);
#     2 the mean r with the four H2 neighbours (Observing, SBI, BFI-2 Aesthetic
#       Sensitivity, Brief MAIA-2 Noticing three-item), as estimated in the H3d model;
#     3 H3d's multiple R.
#   Agreement models: (i) registered, the agreement factor uncorrelated with every
#   validity factor; (ii) the agreement factor free to correlate with them (04's ri_free).
#   Benchmarks: the registered 8, the 10, and the 7 without U4, through the same code.
#   Evaluation (search half, test half, full sample) in lavaan: the H1 rule and its
#   quantities, omega, alpha; latent r with every H2 neighbour from its measure's own
#   model as 04 and 15 build them (FFMQ five facets, SBI, BFI-2 three facets, Brief MAIA-2
#   eight subscales with decision 46's bound on maia_12), the H3d model's R; latent r with
#   past-week pause days and the vignette self-rating (single-indicator factors, 04's
#   treatment); H6's d (current vs never, 07's sf2b_kg on the observed score, all but one
#   statement answered). The Pareto front of r with Observing against omega in the primary
#   search half. A bootstrap of the test half for the difference between each chosen
#   subset and the 8.
#
# Settings as registered (Section 5): ordinal items, WLSMV, theta, missing = "pairwise",
# conditional.x = FALSE. CHOICES below are written to the JSON.
#
# Writes output/<tag>/21_*.csv and 21_exploratory_subsets.json (aggregates only).
# Run from the repo root:  SF2B_MODE=real SF2B_CORES=10 Rscript study2b/analysis/21_exploratory_subsets.R

# Set Up ---------
## Load libraries ---------
packages <- c("jsonlite", "data.table", "lavaan", "psych", "parallel")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

## Settings ---------
source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")
stopifnot(SF2B$mode == "real")
T0 <- Sys.time()

# The pipeline's own data: 04_validity.R's top-level expressions up to its first fit (as
# 10-18 do): A, D, M, S, DAT, score_items; nothing is fitted.
ex04 <- parse("study2b/analysis/04_validity.R", keep.source = FALSE)
first_fit <- which(vapply(ex04, function(e) grepl("^sf2b_log\\(\"fitting", deparse(e)[1]), TRUE))
stopifnot(length(first_fit) == 1)
for (i in seq_len(first_fit - 1)) eval(ex04[[i]], envir = globalenv())
stopifnot(exists("DAT"), exists("score_items"), !exists("RES"))
ex10 <- parse("study2b/analysis/10_deviation_maia.R", keep.source = FALSE)
is_def <- function(e, nm) is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(nm))
eval(Filter(function(e) is_def(e, "improper_of"), as.list(ex10))[[1]], envir = globalenv())
improper_list <- function(f) { x <- improper_of(f)$list; x[!grepl("^std loading (L_[[:alnum:]_]+|grp) =~ ", x)] }
stopifnot(identical(M$participant, S$participant), nrow(DAT) == nrow(M), identical(rownames(A$vignettes), M$participant))
DAT$vig_self <- A$vignettes[, "vig_self"]

X <- A$scale_raw
I8 <- D$sf8; REV <- D$reversed; SRC <- D$source
ST <- c(I8, "sf2b_08", "sf2b_12", "sf2b_05", "sf2b_11", "sf2b_13")
stopifnot(setequal(ST, D$items), SRC[["sf2b_07"]] == "U4", SRC[["sf2b_08"]] == "A9", SRC[["sf2b_12"]] == "sfc_24",
          SRC[["sf2b_05"]] == "N6", SRC[["sf2b_11"]] == "sf6_44", SRC[["sf2b_13"]] == "sfc_31")
BENCH <- list(the_8 = I8, the_10 = D$sf10, the_7_without_U4 = setdiff(I8, "sf2b_07"))
FFMQ <- SF2B_MEASURES$ffmq
H2 <- SF2B_H2
VERSIONS <- c(i = FALSE, ii = TRUE)                       # agreement free to correlate with the validity factors?
VERSION_LABEL <- c(i = "registered: agreement uncorrelated with validity factors", ii = "agreement free to correlate with validity factors")
OBJECTIVES <- c(r_obs = "r with FFMQ Observing", mean4 = "mean r with the four H2 neighbours (H3d model)", R = "H3d multiple R")
OMEGA_MIN <- .75
MAIA_BOUND <- 9                                           # decision 46 under theta: l^2 <= 9, std loading <= .949
SEED <- 20261116L
N_SPLITS <- as.integer(Sys.getenv("SF21_SPLITS", "20"))
NEED <- 5L                                                # passers found per objective, version and split
CAP <- 400L                                               # at most this many lavaan checks per objective, version and split
B_BOOT <- as.integer(Sys.getenv("SF21_BOOT", "300"))
MATERIAL <- .05                                           # a material separation: r with Observing at least .05 below the 8's
CONV_DROP_R <- .05; CONV_DROP_D <- .10                    # convergent-validity flags (against the 8, same sample)
CORES <- SF2B$cores
g_cur <- fifelse(M$practice_group == "current", 1L, fifelse(M$practice_group == "never", 0L, NA_integer_))
key <- function(it) ifelse(it %in% REV, -1, 1)
lab <- function(items) paste(SRC[ST[ST %in% items]], collapse = " ")
order_items <- function(items) items[order(items %in% REV)]   # S takes its sign from an unreversed statement

## Candidates ---------
SUBSETS <- list()
for (k in 5:10) for (cmb in combn(13, k, simplify = FALSE)) {
  it <- ST[cmb]
  if (any(it %in% REV) && any(!it %in% REV)) SUBSETS[[length(SUBSETS) + 1L]] <- it
}
SUB_LABEL <- vapply(SUBSETS, lab, "")
stopifnot(length(SUBSETS) == 6971)
# SF21_MAX_SUBSETS (timing runs only): a random sample of candidates plus the benchmarks,
# outputs written as 21dry_*. Unset, every candidate.
MAXS <- as.integer(Sys.getenv("SF21_MAX_SUBSETS", "0"))
if (MAXS > 0) { set.seed(1); keep <- sort(unique(c(match(vapply(BENCH, lab, ""), SUB_LABEL), sample(length(SUBSETS), MAXS))))
  SUBSETS <- SUBSETS[keep]; SUB_LABEL <- SUB_LABEL[keep] }
out21 <- function(f) sf2b_out(if (MAXS > 0) sub("^21_", "21dry_", f) else f)
SUB_ID <- seq_along(SUBSETS)
bench_id <- vapply(BENCH, function(b) which(SUB_LABEL == lab(b)), 0L)
sf2b_log("21: ", length(SUBSETS), " candidate subsets (5-10 statements, both keying directions)")

## Splits ---------
# Stratified by practice group (current, never, past or occasional, missing): within each
# stratum floor(n/2) to the search half, the odd strata's spare places assigned at random
# so the halves are 360 and 360.
strata <- fifelse(is.na(M$practice_group), "missing", M$practice_group)
set.seed(SEED)
SPLITS <- lapply(seq_len(N_SPLITS), function(s) {
  ns <- table(strata); base <- floor(ns / 2); spare <- nrow(M) / 2 - sum(base)
  odd <- names(ns)[ns %% 2 == 1]
  extra <- if (spare > 0) sample(odd, spare) else character()
  search <- sort(unlist(lapply(names(ns), function(g) {
    ix <- which(strata == g); n <- base[[g]] + (g %in% extra); ix[sample.int(length(ix), n)] })))
  list(search = search, test = setdiff(seq_len(nrow(M)), search))
})
stopifnot(all(vapply(SPLITS, function(x) length(x$search) == 360 && length(x$test) == 360, TRUE)))
FULL <- seq_len(nrow(M))

CHOICES <- list(
  status = "EXPLORATORY, NOT REGISTERED, POST HOC: specified after the confirmatory run and analyses 11-19 were seen. A subset chosen here is a hypothesis for a new sample.",
  pool = "The 13 statements: the 8 (N4, sf6_36 R, V10, sf6_27, sfc_36, sfc_x2 R, sfc_35, U4 R), the 10's A9 and sfc_24, and the unscored reversed N6, sfc_31, sf6_44. Repository ids only.",
  candidates = "Every subset of 5 to 10 statements with at least one reversed and at least one unreversed statement: 6,971 (7,007 subsets of 5-10, less 29 with no reversed and 7 with no unreversed statement). The brief asked for at least one reversed; the 7 all-reversed subsets are also left out because the agreement factor (loadings 1 on the unrecoded statements) cannot be told from content when every statement runs the same way.",
  split = sprintf("Halves of 360, stratified by practice group (current, never, past or occasional, missing); floor(n/2) of each stratum to the search half, the odd strata's spare places at random. %d splits drawn in order after set.seed(%d); split 1 is the primary split.", N_SPLITS, SEED),
  screen = "Point estimates for every candidate from the search half's pairwise polychoric matrix (lavCor, missing = pairwise) and lavaan's DWLS weights for it, by a direct DWLS fit (nlminb, analytic gradient) of the same theta-parameterized models: the H1 model, 04's FFMQ model (five correlated facets) and 04's H3d model (the four H2 neighbours, correlated). Thresholds are saturated and the weight matrix diagonal, so the estimates equal lavaan's WLSMV point estimates; the agreement factor's variance is kept positive and correlations in (-1, 1), so a candidate whose lavaan solution would be improper can differ, which the lavaan check catches. Agreement with lavaan is reported (21_screen_check.csv).",
  constraints = "Checked in lavaan on the search half: the H1 rule exactly as registered (sf2b_h1 with the agreement factor as registered); content omega >= .75 (03's formula); no improper solution (10's improper_of(): std loading >= 1, negative variances, factor correlation >= 1, non-positive-definite factor covariance) in the H1 model, the FFMQ model or the H3d model of that agreement version; every factor oriented as keyed (sf2b_orientation).",
  selection = sprintf("Per split, agreement version and objective: the candidates with screened omega >= %.2f are walked in order of the screened objective, each checked in lavaan, until the %d best meeting every constraint are found (at most %d checked). The best is the chosen subset; the five are its runners-up. Statements are entered unreversed first so S takes its sign from an unreversed statement.", OMEGA_MIN, NEED, CAP),
  objectives = "Minimised: (1) latent r of S with FFMQ Observing in 04's FFMQ model; (2) the mean of S's latent r with the four H2 neighbours (Observing, SBI, BFI-2 Aesthetic Sensitivity, Brief MAIA-2 Noticing) in the H3d model; (3) the H3d model's multiple R = sqrt(r' Sigma^-1 r). Objective 2 is read from the H3d model, which holds all four (Noticing without maia_12), so the search needs no bound; the evaluation tables also give the four cells from each measure's own model (04/15's: Brief MAIA-2 with decision 46's bound).",
  agreement_versions = "(i) registered: the agreement factor uncorrelated with every validity factor; (ii) the agreement factor free to correlate with every validity factor in the model (04's ri_free), as in 14's sensitivity analysis. The H1 model, omega and the H1 rule do not involve validity factors, so they are the same in both.",
  evaluation = "Search half, test half and full sample, in lavaan with SEs: sf2b_h1, 03's alpha; 04/15's measure models for the four H2 neighbours; the H3d model's R; single-indicator latent r with past-week pause days and the vignette self-rating (04's treatment; the agreement factor uncorrelated or free per version); H6's d, current vs never on the observed mean (07's sf2b_kg, Student's t, pooled SD), all but one statement answered (the registered rule's form: 7 of 8, 9 of 10).",
  convergent_flags = sprintf("A subset is flagged as gaining separation by losing convergent validity when, in the same sample, its r with pause days or with the vignette self-rating is at least %.2f below the 8's, or its H6 d is at least %.2f below the 8's or below .50 (H6's registered bound).", CONV_DROP_R, CONV_DROP_D),
  material = sprintf("Material separation in the test half: r with Observing at least %.2f below the registered 8's in the same half and agreement version, with the H1 rule met, omega >= .75, no improper solution and no convergent-validity flag in that half. The bootstrap CI of the difference (%d test-half resamples, percentile) is reported beside it.", MATERIAL, B_BOOT),
  pareto = "Primary split, search half: the subsets not dominated on (lower r with Observing, higher omega) among the candidates meeting the H1 rule and with no improper solution, omega unrestricted. Built iteratively: the screened front is checked in lavaan, failures removed, the front recomputed, until every member passes.",
  stability = "The search repeated in each split; for each agreement version and objective the chosen subset and the five best are recorded. Reported: how often each subset is chosen; each statement's share of the chosen subsets and of the five best; the chosen subset's gain over the 8 in its search half and in its own test half (the optimism of selection).",
  bootstrap = sprintf("Primary split's test half, %d resamples (set.seed(%d)): the difference between each chosen subset and the 8 (objective 1 from the FFMQ model; 2 and 3 from the H3d model), both fitted in each resample (se and test off). Failed or improper fits are left out and counted.", B_BOOT, SEED + 1L),
  settings = "Ordinal, WLSMV, theta, missing = pairwise, conditional.x = FALSE (sf2b_sem). Aggregates only; no participant identifiers.")

# DWLS screen ---------
# Direct DWLS on a polychoric matrix, theta parameterization: Sigma* = L Phi L' + I, the
# model correlation matrix its standardization, loss sum_{i<j} w_ij (r_ij - rho_ij)^2.
# Factors: S (variance 1), RI (loadings 1 on the statements, variance psi), then the
# validity factors (variance 1, correlated with S and each other; with RI if ri_free).
dwls_fit <- function(R, W, items, vfac = list(), ri_free = FALSE, vstart = NULL) {
  vv <- unlist(vfac, use.names = FALSE); vars <- c(items, vv)
  k <- length(items); m <- length(vv); nF <- length(vfac); nf <- 2L + nF
  Rs <- R[vars, vars]; Ws <- W[vars, vars]; lt <- lower.tri(Rs)
  fidx <- rep(seq_len(nF), lengths(vfac))
  lidx <- rbind(cbind(seq_len(k), 1L), if (m) cbind(k + seq_len(m), 2L + fidx))
  nl <- nrow(lidx)
  ff <- if (nF > 1) t(combn(nF, 2)) else matrix(0L, 0, 2); nff <- nrow(ff)
  ffa <- cbind(2L + ff[, 1], 2L + ff[, 2]); ffb <- ffa[, 2:1, drop = FALSE]
  L0 <- matrix(0, k + m, nf); L0[seq_len(k), 2] <- 1
  o_psi <- nl + 1L; o_sf <- o_psi + seq_len(nF); o_ff <- o_psi + nF + seq_len(nff)
  o_rf <- o_psi + nF + nff + seq_len(if (ri_free) nF else 0L)
  vf <- if (nF) 3:nf else integer()
  build <- function(p) {
    Lm <- L0; Lm[lidx] <- p[seq_len(nl)]
    psi <- exp(p[o_psi]); Ph <- diag(nf); Ph[2, 2] <- psi
    if (nF) { r <- tanh(p[o_sf]); Ph[1, vf] <- r; Ph[vf, 1] <- r }
    if (nff) { r <- tanh(p[o_ff]); Ph[ffa] <- r; Ph[ffb] <- r }
    if (ri_free && nF) { r <- tanh(p[o_rf]) * sqrt(psi); Ph[2, vf] <- r; Ph[vf, 2] <- r }
    list(L = Lm, Ph = Ph, psi = psi)
  }
  implied <- function(b) { Sg <- tcrossprod(b$L %*% b$Ph, b$L); diag(Sg) <- diag(Sg) + 1; d <- 1 / sqrt(diag(Sg)); list(Sg = Sg, d = d, Rh = Sg * outer(d, d)) }
  obj <- function(p) { im <- implied(build(p)); e <- (Rs - im$Rh)[lt]; sum(Ws[lt] * e * e) }
  grad <- function(p) {
    b <- build(p); im <- implied(b); d <- im$d
    E <- -2 * Ws * (Rs - im$Rh); diag(E) <- 0
    H <- 0.5 * E * outer(d, d); diag(H) <- -0.5 * rowSums(E * im$Rh) * d^2
    LP <- b$L %*% b$Ph
    GL <- 2 * H %*% LP; GP <- crossprod(b$L, H %*% b$L)
    g <- numeric(o_psi + nF + nff + length(o_rf)); g[seq_len(nl)] <- GL[lidx]
    gpsi <- GP[2, 2] * b$psi
    if (ri_free && nF) gpsi <- gpsi + sum(GP[2, vf] * b$Ph[2, vf])
    g[o_psi] <- gpsi
    if (nF) g[o_sf] <- 2 * GP[1, vf] * (1 - b$Ph[1, vf]^2)
    if (nff) g[o_ff] <- 2 * GP[ffa] * (1 - b$Ph[ffa]^2)
    if (ri_free && nF) { r <- b$Ph[2, vf] / sqrt(b$psi); g[o_rf] <- 2 * GP[2, vf] * (1 - r^2) * sqrt(b$psi) }
    g
  }
  kv <- c(key(items), if (m) ifelse(vv %in% REV_VAL, -1, 1))
  st <- c(kv[lidx[, 1]], log(.3), rep(atanh(.6), nF), rep(atanh(.4), nff), rep(0, length(o_rf)))
  if (!is.null(vstart)) { j <- k + seq_len(m); st[j] <- vstart$lv[vv]; if (nff) st[o_ff] <- vstart$ff }
  o <- nlminb(st, obj, grad, control = list(iter.max = 1000, eval.max = 2000, rel.tol = 1e-11))
  b <- build(o$par); im <- implied(b)
  sdv <- 1 / im$d; Ls <- b$L %*% diag(sqrt(diag(b$Ph))) / sdv; dimnames(Ls) <- list(vars, c("S", "RI", names(vfac)))
  Ph <- b$Ph; dimnames(Ph) <- list(colnames(Ls), colnames(Ls)); Cr <- cov2cor(Ph)
  gmax <- max(abs(grad(o$par)))
  ev <- min(eigen(Ph, symmetric = TRUE, only.values = TRUE)$values)
  oriented <- sum(Ls[items, "S"] * key(items)) > 0 &&
    all(vapply(names(vfac), function(f) sum(Ls[vfac[[f]], f] * ifelse(vfac[[f]] %in% REV_VAL, -1, 1)) > 0, TRUE))
  res <- (Rs - im$Rh)[items, items]
  list(ok = gmax < 1e-3 && oriented, gmax = gmax, oriented = oriented, boundary = b$psi < 1e-3 || any(abs(Cr[upper.tri(Cr)]) > .995) || ev <= 1e-8,
       cor = Cr, Lstd = Ls, theta_std = im$d^2, maxres = max(abs(res[lower.tri(res)])), par = o$par, lv = o$par[k + seq_len(m)],
       ff = if (nff) o$par[o_ff] else NULL)
}
REV_VAL <- unlist(lapply(D$validity, `[[`, "reverse"))
dwls_omega <- function(f, items) { k <- key(items); l <- f$Lstd[items, "S"] * k; r <- f$Lstd[items, "RI"] * k
  sum(l)^2 / (sum(l)^2 + sum(r)^2 + sum(f$theta_std[seq_along(items)])) }
mult_R <- function(C) { r <- C["S", H2]; Sig <- C[H2, H2]; v <- tryCatch(drop(t(r) %*% solve(Sig, r)), error = function(e) NA_real_); if (is.na(v) || v < 0) NA_real_ else sqrt(v) }
VF_FFMQ <- setNames(lapply(FFMQ, function(s) score_items[[s]]), FFMQ)
VF_H2 <- setNames(lapply(H2, function(s) score_items[[s]]), H2)
PC_VARS <- unique(c(ST, unlist(VF_FFMQ), unlist(VF_H2)))

polychoric_of <- function(rows) {
  lc <- lavCor(DAT[rows, PC_VARS], ordered = PC_VARS, missing = "pairwise", output = "fit", se = "robust.sem", estimator = "DWLS", do.fit = FALSE)
  R <- lavInspect(lc, "sampstat")$cov; wobs <- lavInspect(lc, "wls.obs"); Wd <- diag(lavInspect(lc, "wls.v"))
  nv <- nrow(R); nth <- length(wobs) - nv * (nv - 1) / 2; ix <- which(lower.tri(R), arr.ind = TRUE)
  stopifnot(max(abs(wobs[(nth + 1):length(wobs)] - R[ix])) < 1e-12)
  W <- matrix(0, nv, nv, dimnames = dimnames(R)); W[ix] <- Wd[(nth + 1):length(wobs)]; W <- W + t(W)
  list(R = R, W = W)
}
screen_all <- function(PC) {
  ref5 <- lapply(VERSIONS, function(v) dwls_fit(PC$R, PC$W, I8, VF_FFMQ, v))
  ref3 <- lapply(VERSIONS, function(v) dwls_fit(PC$R, PC$W, I8, VF_H2, v))
  rbindlist(lapply(SUB_ID, function(id) {
    it <- order_items(SUBSETS[[id]])
    h <- dwls_fit(PC$R, PC$W, it)
    row <- list(id = id, k = length(it), omega = dwls_omega(h, it), maxres = h$maxres, h1_ok = h$ok, h1_boundary = h$boundary)
    for (vn in names(VERSIONS)) {
      v <- VERSIONS[[vn]]
      f5 <- dwls_fit(PC$R, PC$W, it, VF_FFMQ, v, vstart = list(lv = setNames(ref5[[vn]]$lv, unlist(VF_FFMQ)), ff = ref5[[vn]]$ff))
      f3 <- dwls_fit(PC$R, PC$W, it, VF_H2, v, vstart = list(lv = setNames(ref3[[vn]]$lv, unlist(VF_H2)), ff = ref3[[vn]]$ff))
      row[[paste0("r_obs_", vn)]] <- f5$cor["S", "ffmq_observing"]
      row[[paste0("mean4_", vn)]] <- mean(f3$cor["S", H2])
      row[[paste0("R_", vn)]] <- mult_R(f3$cor)
      row[[paste0("ok_", vn)]] <- f5$ok && f3$ok
      row[[paste0("boundary_", vn)]] <- f5$boundary || f3$boundary
    }
    row
  }))
}

# lavaan models ---------
alpha_of <- function(items, rows) { Z <- X[rows, items, drop = FALSE]; rv <- intersect(items, REV); Z[, rv] <- 7 - Z[, rv]   # 03's
  suppressWarnings(psych::alpha(Z[complete.cases(Z), ], check.keys = FALSE, warnings = FALSE)$total$raw_alpha) }
check_h1 <- function(items, rows) {
  it <- order_items(items)
  h <- sf2b_try(sf2b_h1(X[rows, , drop = FALSE], it, REV))
  if (inherits(h, "error") || !isTRUE(h$converged))
    return(list(converged = FALSE, error = if (inherits(h, "error")) conditionMessage(h) else h$error, supported = FALSE, pass = FALSE))
  imp <- paste(improper_list(h$fit), collapse = "; ")
  bad <- sf2b_orientation(h$fit, it, character(), D)
  list(converged = TRUE, error = if (length(bad)) "S not oriented as keyed" else "", supported = h$supported, cfi = h$stats$cfi.scaled,
       rmsea = h$stats$rmsea.scaled, srmr = h$stats$srmr, max_residual = h$max_residual,
       max_residual_pair = paste(SRC[strsplit(h$max_residual_pair, " ~~ ")[[1]]], collapse = " ~~ "),
       cfi_ok = h$cfi_ok, local_ok = h$local_ok, close_ok = h$close_ok, closefit_p = h$closefit_p, omega = h$omega,
       agreement_share = h$ri_share, improper = imp,
       pass = isTRUE(h$supported) && h$omega >= OMEGA_MIN && !nzchar(imp) && !length(bad))
}
SINGLE1 <- c("crit_pause_days", "vig_self")
BLOCKS <- list(ffmq = FFMQ, sbi = "sbi_moment", bfi = SF2B_MEASURES$bfi, maia = SF2B_MEASURES$maia, h3d = H2,
               crit_pause_days = "crit_pause_days", vig_self = "vig_self")
val_model <- function(items, b, ri_free) {
  syn <- sf2b_scale_syntax(items); ord <- items; scores <- character()
  if (b %in% SINGLE1) {
    syn <- c(syn, paste0("L_", b, " =~ 1*", b), paste0(b, " ~~ 0*", b), paste0("L_", b, " ~~ 1*L_", b),
             paste0("S ~~ r_", b, "*L_", b), paste0("RI ~~ ", if (ri_free) paste0("ri_", b, "*") else "0*", "L_", b))
    ord <- c(ord, b)
  } else {
    scores <- BLOCKS[[b]]
    syn <- c(syn, sf2b_validity_syntax(scores, D, ri_free = ri_free))
    if ("maia_notworrying" %in% scores) syn <- c(syn, "maia_notworrying =~ l_maia12*maia_12", sprintf("l_maia12^2 < %s", MAIA_BOUND))
    ord <- c(ord, unlist(score_items[scores]))
  }
  list(syn = syn, ord = unique(ord), scores = scores)
}
fit_val <- function(items, b, ri_free, rows, fast = FALSE) {
  it <- order_items(items); m <- val_model(it, b, ri_free)
  f <- sf2b_try(sf2b_sem(m$syn, DAT[rows, , drop = FALSE], ordered = m$ord, fast = fast))
  if (sf2b_failed(f)) return(list(ok = FALSE, error = sf2b_msg(f), improper = ""))
  bad <- sf2b_orientation(f, it, m$scores, D)
  imp <- paste(improper_list(f), collapse = "; ")
  C <- lavInspect(f, "cor.lv")
  out <- list(ok = !length(bad), error = if (length(bad)) paste("not oriented as keyed:", paste(bad, collapse = ", ")) else "", improper = imp,
              cor_S = C["S", setdiff(colnames(C), c("S", "RI"))])
  names(out$cor_S) <- sub("^L_", "", names(out$cor_S))
  if (b == "h3d") out$R <- mult_R(C)
  if (!fast) {
    pe <- parameterEstimates(f); r <- pe[grepl("^r_", pe$label), ]
    out$ci <- data.table(score = sub("^r_", "", r$label), r = r$est, lower = r$ci.lower, upper = r$ci.upper)
    fs <- sf2b_fit_stats(f); out$fit <- fs[c("cfi.scaled", "rmsea.scaled", "srmr")]
  }
  out
}
check_val <- function(items, rows, vn) {
  v <- VERSIONS[[vn]]
  a <- fit_val(items, "ffmq", v, rows, fast = TRUE); b <- fit_val(items, "h3d", v, rows, fast = TRUE)
  pass <- a$ok && b$ok && !nzchar(a$improper) && !nzchar(b$improper)
  list(pass = pass, r_obs = if (a$ok) a$cor_S[["ffmq_observing"]] else NA_real_,
       mean4 = if (b$ok) mean(b$cor_S[H2]) else NA_real_, R = if (b$ok) b$R else NA_real_,
       why = paste(c(if (!a$ok) paste("FFMQ:", a$error), if (nzchar(a$improper)) paste("FFMQ improper:", a$improper),
                     if (!b$ok) paste("H3d:", b$error), if (nzchar(b$improper)) paste("H3d improper:", b$improper)), collapse = "; "))
}

# Search in one split ---------
search_split <- function(s) {
  t0 <- Sys.time(); sp <- SPLITS[[s]]
  PC <- polychoric_of(sp$search)
  scr <- screen_all(PC)
  t_screen <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  H1C <- new.env(); VC <- new.env(); CHK <- list()
  h1_of <- function(id) { k <- as.character(id); if (is.null(H1C[[k]])) H1C[[k]] <- check_h1(SUBSETS[[id]], sp$search); H1C[[k]] }
  val_of <- function(id, vn) { k <- paste(id, vn); if (is.null(VC[[k]])) VC[[k]] <- check_val(SUBSETS[[id]], sp$search, vn); VC[[k]] }
  chosen <- list(); walk <- list()
  for (vn in names(VERSIONS)) for (ob in names(OBJECTIVES)) {
    col <- paste0(ob, "_", vn)
    cand <- scr[omega >= OMEGA_MIN & h1_ok == TRUE & get(paste0("ok_", vn)) == TRUE][order(get(col))]
    found <- 0L; n <- 0L; fails <- c(h1_rule = 0L, omega = 0L, improper_h1 = 0L, validity = 0L)
    for (id in cand$id) {
      if (found >= NEED || n >= CAP) break
      n <- n + 1L; h <- h1_of(id)
      if (!isTRUE(h$supported)) { fails["h1_rule"] <- fails["h1_rule"] + 1L; next }
      if (!isTRUE(h$omega >= OMEGA_MIN)) { fails["omega"] <- fails["omega"] + 1L; next }
      if (!isTRUE(h$pass)) { fails["improper_h1"] <- fails["improper_h1"] + 1L; next }
      vv <- val_of(id, vn)
      CHK[[length(CHK) + 1L]] <- data.table(split = s, version = vn, id = id, r_obs_dwls = scr[[paste0("r_obs_", vn)]][id],
                                            r_obs_lavaan = vv$r_obs, mean4_dwls = scr[[paste0("mean4_", vn)]][id], mean4_lavaan = vv$mean4,
                                            R_dwls = scr[[paste0("R_", vn)]][id], R_lavaan = vv$R,
                                            omega_dwls = scr$omega[id], omega_lavaan = h$omega)
      if (!vv$pass) { fails["validity"] <- fails["validity"] + 1L; next }
      found <- found + 1L
      chosen[[length(chosen) + 1L]] <- data.table(split = s, version = vn, objective = ob, rank = found, id = id, k = length(SUBSETS[[id]]),
                                                  statements = SUB_LABEL[id], value_search = vv[[ob]], r_obs = vv$r_obs, mean4 = vv$mean4, R = vv$R,
                                                  omega = h$omega, cfi = h$cfi, max_residual = h$max_residual)
    }
    walk[[length(walk) + 1L]] <- data.table(split = s, version = vn, objective = ob, eligible_screened = nrow(cand), checked = n, found = found,
                                            fail_h1_rule = fails[["h1_rule"]], fail_omega = fails[["omega"]], fail_improper_or_orientation_h1 = fails[["improper_h1"]],
                                            fail_validity_model = fails[["validity"]])
  }
  chosen <- rbindlist(chosen); walk <- rbindlist(walk)
  # Each chosen subset and the 8 in this split's test half (and the 8 in the search half).
  tst <- list(); EC <- new.env()
  ev <- function(id, half, vn, ob, role) {
    k <- paste(id, half, vn)
    if (is.null(EC[[k]])) { rows <- sp[[half]]; h1k <- paste(id, half)
      if (is.null(EC[[h1k]])) EC[[h1k]] <- check_h1(SUBSETS[[id]], rows)
      EC[[k]] <- list(h = EC[[h1k]], v = check_val(SUBSETS[[id]], rows, vn)) }
    h <- EC[[k]]$h; vv <- EC[[k]]$v
    data.table(split = s, half = half, version = vn, objective = ob, role = role, id = id, h1_supported = isTRUE(h$supported), omega = h$omega %||% NA_real_,
               proper = isTRUE(h$pass) && vv$pass, value = vv[[ob]])
  }
  for (vn in names(VERSIONS)) for (ob in names(OBJECTIVES)) {
    best <- chosen[version == vn & objective == ob & rank == 1]
    tst[[length(tst) + 1L]] <- ev(bench_id[["the_8"]], "search", vn, ob, "the 8")
    tst[[length(tst) + 1L]] <- ev(bench_id[["the_8"]], "test", vn, ob, "the 8")
    if (nrow(best)) tst[[length(tst) + 1L]] <- ev(best$id, "test", vn, ob, "chosen")
  }
  list(split = s, screen = if (s == 1L) scr else scr[, .(id, omega, r_obs_i, mean4_i, R_i, r_obs_ii, mean4_ii, R_ii)], chosen = chosen, walk = walk,
       test = rbindlist(tst), check = rbindlist(CHK), polychoric = if (s == 1L) PC else NULL,
       h1_checked = as.list(H1C), secs = c(screen = t_screen, total = as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

# Cluster ---------
cl <- makeCluster(CORES)
invisible(clusterCall(cl, function(wd, env) { setwd(wd); do.call(Sys.setenv, as.list(env)); NULL }, getwd(),
            c(SF2B_MODE = SF2B$mode, SF2B_DATA_DIR = SF2B$data_dir, SF2B_CORES = "1")))
invisible(clusterEvalQ(cl, { suppressMessages({ library(data.table); library(lavaan); library(psych)
  source("study2b/analysis/config.R"); source("study2b/analysis/sf2b_models.R") }); NULL }))
EXPORT <- c("DAT", "X", "D", "score_items", "I8", "REV", "SRC", "ST", "SUBSETS", "SUB_LABEL", "SUB_ID", "bench_id", "SPLITS", "FFMQ", "H2",
            "VERSIONS", "OBJECTIVES", "OMEGA_MIN", "MAIA_BOUND", "NEED", "CAP", "key", "lab", "order_items", "dwls_fit", "REV_VAL", "dwls_omega",
            "mult_R", "VF_FFMQ", "VF_H2", "PC_VARS", "polychoric_of", "screen_all", "alpha_of", "check_h1", "SINGLE1", "BLOCKS", "val_model",
            "fit_val", "check_val", "search_split", "improper_of", "improper_list", "FULL", "g_cur")
clusterExport(cl, EXPORT)

# Reference check of the screen against lavaan (full sample, the 8) ---------
PCF <- polychoric_of(FULL)
ref_chk <- rbindlist(lapply(names(VERSIONS), function(vn) {
  f5 <- dwls_fit(PCF$R, PCF$W, order_items(I8), VF_FFMQ, VERSIONS[[vn]]); f3 <- dwls_fit(PCF$R, PCF$W, order_items(I8), VF_H2, VERSIONS[[vn]])
  lv <- check_val(I8, FULL, vn)
  data.table(version = vn, r_obs_dwls = f5$cor["S", "ffmq_observing"], r_obs_lavaan = lv$r_obs, mean4_dwls = mean(f3$cor["S", H2]),
             mean4_lavaan = lv$mean4, R_dwls = mult_R(f3$cor), R_lavaan = lv$R)
}))
h8 <- dwls_fit(PCF$R, PCF$W, order_items(I8)); ref_chk[, `:=`(omega_dwls = dwls_omega(h8, order_items(I8)), omega_lavaan = check_h1(I8, FULL)$omega)]
print(ref_chk)
stopifnot(max(abs(c(ref_chk$r_obs_dwls - ref_chk$r_obs_lavaan, ref_chk$R_dwls - ref_chk$R_lavaan, ref_chk$omega_dwls - ref_chk$omega_lavaan))) < 1e-3)

# Phase 1: the search in every split ---------
sf2b_log("21: phase 1, ", N_SPLITS, " splits on ", CORES, " workers")
t1 <- Sys.time()
SPL <- parLapplyLB(cl, seq_len(N_SPLITS), search_split)
secs_search <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
sf2b_log("21: phase 1 done in ", round(secs_search), " s")
CHOSEN <- rbindlist(lapply(SPL, `[[`, "chosen")); WALK <- rbindlist(lapply(SPL, `[[`, "walk")); TEST <- rbindlist(lapply(SPL, `[[`, "test"))
CHECK <- rbindlist(lapply(SPL, `[[`, "check"))
SCR1 <- SPL[[1]]$screen
SPLIT_SECS <- rbindlist(lapply(SPL, function(x) data.table(split = x$split, screen_secs = x$secs[["screen"]], total_secs = x$secs[["total"]])))

# Phase 2: the Pareto front (primary split, search half) ---------
pareto <- function(dt, rcol) { d <- dt[order(-omega, get(rcol))]; best <- Inf; keep <- logical(nrow(d))
  for (i in seq_len(nrow(d))) { if (d[[rcol]][i] < best) { keep[i] <- TRUE; best <- d[[rcol]][i] } }; d[keep] }
H1P <- SPL[[1]]$h1_checked                                  # lavaan H1 checks already made in split 1
P1 <- SPLITS[[1]]$search
front <- list()
for (vn in names(VERSIONS)) {
  rcol <- paste0("r_obs_", vn); excl <- integer(); VP <- list()
  repeat {
    pool <- SCR1[h1_ok == TRUE & get(paste0("ok_", vn)) == TRUE & !(id %in% excl)]
    fr <- pareto(pool, rcol)
    need_h1 <- setdiff(fr$id, as.integer(names(H1P)))
    if (length(need_h1)) { r <- parLapplyLB(cl, need_h1, function(id) check_h1(SUBSETS[[id]], SPLITS[[1]]$search)); names(r) <- need_h1; H1P[names(r)] <- r }
    h1_fail <- fr$id[!vapply(as.character(fr$id), function(k) isTRUE(H1P[[k]]$supported) && !nzchar(H1P[[k]]$improper %||% "x") && !nzchar(H1P[[k]]$error %||% ""), TRUE)]
    need_v <- setdiff(setdiff(fr$id, h1_fail), as.integer(names(VP)))
    if (length(need_v)) { r <- parLapplyLB(cl, need_v, function(id, vn) check_val(SUBSETS[[id]], SPLITS[[1]]$search, vn), vn = vn); names(r) <- need_v; VP[names(r)] <- r }
    v_fail <- setdiff(fr$id, h1_fail)[!vapply(as.character(setdiff(fr$id, h1_fail)), function(k) isTRUE(VP[[k]]$pass), TRUE)]
    if (!length(h1_fail) && !length(v_fail)) break
    excl <- c(excl, h1_fail, v_fail)
  }
  front[[vn]] <- cbind(version = vn, fr[, .(id, k, omega, r_obs = get(rcol))], rejected_on_the_way = length(excl))
}
FRONT <- rbindlist(front)
FRONT[, statements := SUB_LABEL[id]]
sf2b_log("21: Pareto fronts: ", paste(names(front), vapply(front, nrow, 0L), sep = " = ", collapse = ", "))

# Phase 3: evaluation of the reported subsets ---------
TOPN <- CHOSEN[split == 1]
REPORT_IDS <- unique(c(bench_id, TOPN$id, FRONT$id))
ROLE <- rbindlist(list(data.table(id = bench_id, role = paste("benchmark:", names(BENCH))),
                       TOPN[, .(id, role = sprintf("best %d on %s, version %s", rank, objective, version))],
                       FRONT[, .(id, role = sprintf("Pareto front, version %s", version))]))
ROLE <- ROLE[, .(role = paste(role, collapse = "; ")), by = id]
SAMPLES <- list(search = SPLITS[[1]]$search, test = SPLITS[[1]]$test, full = FULL)
ETASKS <- list()
for (id in REPORT_IDS) for (sm in names(SAMPLES)) {
  ETASKS[[length(ETASKS) + 1L]] <- list(id = id, sample = sm, kind = "h1")
  for (vn in names(VERSIONS)) for (b in names(BLOCKS)) ETASKS[[length(ETASKS) + 1L]] <- list(id = id, sample = sm, kind = b, version = vn)
}
eval_task <- function(t) {
  rows <- SAMPLES[[t$sample]]; t0 <- Sys.time()
  out <- if (t$kind == "h1") { h <- check_h1(SUBSETS[[t$id]], rows); h$alpha <- alpha_of(SUBSETS[[t$id]], rows); h }
         else fit_val(SUBSETS[[t$id]], t$kind, VERSIONS[[t$version]], rows)
  c(list(task = t, secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))), out)
}
clusterExport(cl, c("SAMPLES", "eval_task"))
sf2b_log("21: phase 3, ", length(REPORT_IDS), " subsets, ", length(ETASKS), " models")
t3 <- Sys.time()
EV <- parLapplyLB(cl, ETASKS, eval_task)
secs_eval <- as.numeric(difftime(Sys.time(), t3, units = "secs"))

# Phase 4: bootstrap of the test half ---------
BEST1 <- TOPN[rank == 1]
set.seed(SEED + 1L)
TEST1 <- SPLITS[[1]]$test
BIDX <- lapply(seq_len(B_BOOT), function(b) TEST1[sample.int(length(TEST1), length(TEST1), replace = TRUE)])
boot_task <- function(b) {
  rows <- BIDX[[b]]; out <- list()
  for (vn in names(VERSIONS)) {
    v <- VERSIONS[[vn]]
    ids5 <- unique(c(bench_id[["the_8"]], BEST1[version == vn & objective == "r_obs", id]))
    ids3 <- unique(c(bench_id[["the_8"]], BEST1[version == vn & objective %in% c("mean4", "R"), id]))
    for (id in ids5) { f <- fit_val(SUBSETS[[id]], "ffmq", v, rows, fast = TRUE)
      out[[length(out) + 1L]] <- data.table(b = b, version = vn, id = id, model = "ffmq", ok = f$ok && !nzchar(f$improper),
                                            r_obs = if (f$ok) f$cor_S[["ffmq_observing"]] else NA_real_, mean4 = NA_real_, R = NA_real_) }
    for (id in ids3) { f <- fit_val(SUBSETS[[id]], "h3d", v, rows, fast = TRUE)
      out[[length(out) + 1L]] <- data.table(b = b, version = vn, id = id, model = "h3d", ok = f$ok && !nzchar(f$improper),
                                            r_obs = NA_real_, mean4 = if (f$ok) mean(f$cor_S[H2]) else NA_real_, R = if (f$ok) f$R else NA_real_) }
  }
  rbindlist(out)
}
clusterExport(cl, c("BIDX", "BEST1", "boot_task"))
sf2b_log("21: phase 4, bootstrap ", B_BOOT, " resamples")
t4 <- Sys.time()
BOOT <- rbindlist(parLapplyLB(cl, seq_len(B_BOOT), boot_task))
secs_boot <- as.numeric(difftime(Sys.time(), t4, units = "secs"))
stopCluster(cl)

# Tables ---------
num <- function(DT, k = 3) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; cols <- setdiff(cols, c("id", "k", "split", "rank", "b"))
  DT[, (cols) := lapply(.SD, function(x) round(x, k)), .SDcols = cols]; DT }
H1E <- rbindlist(lapply(Filter(function(r) r$task$kind == "h1", EV), function(r)
  data.table(id = r$task$id, sample = r$task$sample, converged = r$converged, h1_rule_met = isTRUE(r$supported), cfi = r$cfi %||% NA_real_,
             rmsea = r$rmsea %||% NA_real_, srmr = r$srmr %||% NA_real_, max_residual = r$max_residual %||% NA_real_,
             max_residual_pair = r$max_residual_pair %||% NA_character_, closefit_p = r$closefit_p %||% NA_real_,
             omega = r$omega %||% NA_real_, alpha = r$alpha %||% NA_real_, agreement_share = r$agreement_share %||% NA_real_,
             improper_h1 = paste(c(r$error, r$improper), collapse = " ")[1])), fill = TRUE)
VE <- rbindlist(lapply(Filter(function(r) r$task$kind != "h1", EV), function(r) {
  base <- data.table(id = r$task$id, sample = r$task$sample, version = r$task$version, block = r$task$kind, ok = r$ok, error = r$error,
                     improper = r$improper, secs = r$secs)
  if (!r$ok || is.null(r$ci)) return(base)
  x <- cbind(base[rep(1, nrow(r$ci))], r$ci)
  if (!is.null(r$R)) x <- rbind(x, cbind(base, data.table(score = "R_h3d", r = r$R, lower = NA_real_, upper = NA_real_)))
  x
}), fill = TRUE)
fwrite(num(copy(VE)), out21("21_model_estimates.csv"))
# One row per subset x sample x version.
cellv <- function(blk, sc) VE[block == blk & score == sc, .(id, sample, version, v = r, lo = lower, hi = upper, bad = !ok | nzchar(improper))]
wide <- Reduce(function(a, b) merge(a, b, by = c("id", "sample", "version"), all = TRUE), list(
  setnames(cellv("ffmq", "ffmq_observing"), c("v", "lo", "hi", "bad"), c("r_obs", "r_obs_lower", "r_obs_upper", "bad_ffmq")),
  setnames(cellv("sbi", "sbi_moment")[, .(id, sample, version, v, bad)], c("v", "bad"), c("r_sbi", "bad_sbi")),
  setnames(cellv("bfi", "bfi_aesthetic")[, .(id, sample, version, v, bad)], c("v", "bad"), c("r_bfi_aesthetic", "bad_bfi")),
  setnames(cellv("maia", "maia_noticing")[, .(id, sample, version, v, bad)], c("v", "bad"), c("r_maia_noticing", "bad_maia")),
  setnames(cellv("h3d", "R_h3d")[, .(id, sample, version, v, bad)], c("v", "bad"), c("R_h3d", "bad_h3d")),
  dcast(VE[block == "h3d" & score %in% H2, .(id, sample, version, score = paste0("h3d_", score), r)], id + sample + version ~ score, value.var = "r"),
  setnames(cellv("crit_pause_days", "crit_pause_days")[, .(id, sample, version, v, bad)], c("v", "bad"), c("r_pause_days", "bad_pause")),
  setnames(cellv("vig_self", "vig_self")[, .(id, sample, version, v, bad)], c("v", "bad"), c("r_vig_self", "bad_vig"))))
wide[, mean4_own_models := (r_obs + r_sbi + r_bfi_aesthetic + r_maia_noticing) / 4]
wide[, mean4_h3d := rowMeans(.SD), .SDcols = paste0("h3d_", H2)]
badc <- grep("^bad_", names(wide), value = TRUE)
for (b in badc) set(wide, which(is.na(wide[[b]])), b, TRUE)   # a missing cell (failed model) counts as not proper
wide[, improper_validity := Reduce(`|`, .SD), .SDcols = badc]
KG <- rbindlist(lapply(REPORT_IDS, function(id) rbindlist(lapply(names(SAMPLES), function(sm) {
  rows <- SAMPLES[[sm]]; it <- SUBSETS[[id]]
  y <- sf2b_score(X[rows, , drop = FALSE], it, REV, need = length(it) - 1L); k <- sf2b_kg(y, g_cur[rows])
  data.table(id = id, sample = sm, h6_d = k$d, h6_d_lower = k$d_lower, h6_d_upper = k$d_upper, n_current = k$n_group, n_never = k$n_reference)
}))))
TAB <- merge(merge(wide, H1E, by = c("id", "sample"), all.x = TRUE), KG, by = c("id", "sample"), all.x = TRUE)
TAB <- merge(TAB, ROLE, by = "id", all.x = TRUE)
TAB[, `:=`(statements = SUB_LABEL[id], k = lengths(SUBSETS)[id])]
ref8 <- TAB[id == bench_id[["the_8"]], .(sample, version, r_obs_8 = r_obs, mean4_8 = mean4_h3d, R_8 = R_h3d, pause_8 = r_pause_days, vig_8 = r_vig_self, d_8 = h6_d, omega_8 = omega)]
TAB <- merge(TAB, ref8, by = c("sample", "version"), all.x = TRUE)
TAB[, `:=`(d_r_obs = r_obs - r_obs_8, d_mean4 = mean4_h3d - mean4_8, d_R = R_h3d - R_8)]
TAB[, convergent_flag := (r_pause_days <= pause_8 - CONV_DROP_R) | (r_vig_self <= vig_8 - CONV_DROP_R) | (h6_d <= d_8 - CONV_DROP_D) | (h6_d < SF2B$h6_d)]
TAB[, sound := h1_rule_met & omega >= OMEGA_MIN & !improper_validity & !nzchar(trimws(improper_h1))]
TAB[, material_separation := sound & !convergent_flag & d_r_obs <= -MATERIAL]
setcolorder(TAB, c("id", "role", "statements", "k", "sample", "version", "h1_rule_met", "omega", "alpha", "cfi", "rmsea", "srmr", "max_residual",
                   "r_obs", "r_obs_lower", "r_obs_upper", "d_r_obs", "mean4_h3d", "d_mean4", "R_h3d", "d_R", "mean4_own_models",
                   "r_sbi", "r_bfi_aesthetic", "r_maia_noticing", "r_pause_days", "r_vig_self", "h6_d", "h6_d_lower", "h6_d_upper",
                   "convergent_flag", "sound", "material_separation"))
TAB <- TAB[order(match(sample, names(SAMPLES)), version, !grepl("benchmark", role), id)]
fwrite(num(copy(TAB)), out21("21_best_subsets.csv"))
FRONT_OUT <- merge(FRONT, TAB[, .(id, sample, version, r_obs_lavaan = r_obs, d_r_obs, omega_lavaan = omega, h1_rule_met, sound, convergent_flag)],
                   by = c("id", "version"), all.x = TRUE, allow.cartesian = TRUE)
fwrite(num(copy(FRONT_OUT[order(version, sample, -omega)])), out21("21_pareto_front.csv"))

# Stability across splits.
fwrite(num(copy(CHOSEN)), out21("21_stability_chosen.csv"))
fwrite(WALK, out21("21_search_walk.csv"))
STAB_SUB <- CHOSEN[rank == 1, .(times_chosen = .N, splits = paste(split, collapse = " "), k = k[1], mean_value_search = mean(value_search)),
                   by = .(version, objective, statements)][order(version, objective, -times_chosen)]
STAB_SUB[, share := times_chosen / N_SPLITS]
fwrite(num(copy(STAB_SUB)), out21("21_stability_subsets.csv"))
STAB_ST <- rbindlist(lapply(c(1L, NEED), function(top) CHOSEN[rank <= top, {
  n_sets <- .N
  data.table(statement = SRC[ST], role = ifelse(ST %in% I8, "the 8", ifelse(ST %in% D$sf10, "the 10", "unscored reversed")),
             reversed = ST %in% REV, share = vapply(ST, function(s) mean(vapply(SUBSETS[id], function(x) s %in% x, TRUE)), 0), sets = n_sets)
}, by = .(version, objective)][, top := top]))
fwrite(num(copy(STAB_ST)), out21("21_stability_statements.csv"))
# The optimism of selection: chosen subset vs the 8, search half and own test half, per split.
TW <- dcast(TEST, split + version + objective ~ role + half, value.var = "value")
setnames(TW, c("chosen_test", "the 8_search", "the 8_test"), c("chosen_test", "the8_search", "the8_test"))
TW <- merge(CHOSEN[rank == 1, .(split, version, objective, statements, k, chosen_search = value_search)], TW, by = c("split", "version", "objective"))
TW <- merge(TW, TEST[role == "chosen", .(split, version, objective, test_h1_met = h1_supported, test_omega = omega, test_proper = proper)], by = c("split", "version", "objective"))
TW[, `:=`(gain_search = chosen_search - the8_search, gain_test = chosen_test - the8_test)]
TW[, material_in_test := objective == "r_obs" & gain_test <= -MATERIAL & test_h1_met & test_omega >= OMEGA_MIN & test_proper]
fwrite(num(copy(TW)), out21("21_stability_splits.csv"))
OPT <- TW[, .(splits = .N, distinct_subsets = uniqueN(statements), mean_gain_search = mean(gain_search), mean_gain_test = mean(gain_test, na.rm = TRUE),
              min_gain_test = min(gain_test, na.rm = TRUE), max_gain_test = max(gain_test, na.rm = TRUE),
              test_h1_met = sum(test_h1_met), test_omega_ok = sum(test_omega >= OMEGA_MIN, na.rm = TRUE),
              gain_test_at_least_material = sum(gain_test <= -MATERIAL, na.rm = TRUE)), by = .(version, objective)]
fwrite(num(copy(OPT)), out21("21_stability_summary.csv"))

# Screen accuracy.
CHECK[, `:=`(d_r_obs = r_obs_lavaan - r_obs_dwls, d_mean4 = mean4_lavaan - mean4_dwls, d_R = R_lavaan - R_dwls, d_omega = omega_lavaan - omega_dwls)]
SCHK <- CHECK[, .(models_compared = .N, max_abs_d_r_obs = max(abs(d_r_obs), na.rm = TRUE), max_abs_d_mean4 = max(abs(d_mean4), na.rm = TRUE),
                  max_abs_d_R = max(abs(d_R), na.rm = TRUE), max_abs_d_omega = max(abs(d_omega), na.rm = TRUE),
                  over_001 = sum(pmax(abs(d_r_obs), abs(d_mean4), abs(d_R), abs(d_omega), na.rm = TRUE) > .001)), by = version]
fwrite(num(copy(SCHK), 5), out21("21_screen_check.csv"))
SCR_OUT <- copy(SCR1)[, statements := SUB_LABEL[id]]
fwrite(num(copy(SCR_OUT), 4), out21("21_screen_primary_search.csv"))

# Bootstrap.
BW <- merge(BOOT[id != bench_id[["the_8"]]], BOOT[id == bench_id[["the_8"]], .(b, version, model, ok8 = ok, r_obs8 = r_obs, mean48 = mean4, R8 = R)], by = c("b", "version", "model"))
BS <- rbindlist(lapply(seq_len(nrow(BEST1)), function(i) {
  x <- BEST1[i]; mdl <- if (x$objective == "r_obs") "ffmq" else "h3d"
  d <- BW[version == x$version & id == x$id & model == mdl]
  dd <- switch(x$objective, r_obs = d$r_obs - d$r_obs8, mean4 = d$mean4 - d$mean48, R = d$R - d$R8)
  okb <- d$ok & d$ok8 & !is.na(dd)
  data.table(version = x$version, objective = x$objective, statements = x$statements, resamples = B_BOOT, used = sum(okb),
             difference_lower = quantile(dd[okb], .025), difference_upper = quantile(dd[okb], .975), boot_mean = mean(dd[okb]))
}))
pt <- TAB[sample == "test", .(id, version, d_r_obs, d_mean4, d_R)]
BS <- merge(BS, merge(BEST1[, .(version, objective, id)], pt, by = c("id", "version")), by = c("version", "objective"))
BS[, difference_test := fifelse(objective == "r_obs", d_r_obs, fifelse(objective == "mean4", d_mean4, d_R))]
BS[, c("d_r_obs", "d_mean4", "d_R") := NULL]
fwrite(num(copy(BS)), out21("21_bootstrap_test_half.csv"))

## Report ---------
rnd <- function(x) if (is.data.frame(x)) x else if (is.list(x)) lapply(x, rnd) else if (is.numeric(x)) signif(x, 4) else x
SMALL <- c("id", "role", "statements", "k", "sample", "version", "h1_rule_met", "omega", "cfi", "rmsea", "srmr", "max_residual", "r_obs", "d_r_obs",
           "mean4_h3d", "R_h3d", "r_pause_days", "r_vig_self", "h6_d", "convergent_flag", "sound", "material_separation")
out <- list(
  analysis = "21: is there a subset of the 13 statements that separates better from FFMQ Observing and the near neighbours while staying a sound scale? EXPLORATORY, NOT REGISTERED, POST HOC (osf.io/d4a5u does not include it). Any subset chosen here is a hypothesis for a new sample, not a result; the registered 8 remains the scale.",
  mode = SF2B$mode, tag = SF2B$tag, written = format(Sys.time(), "%Y-%m-%d %H:%M"),
  software = list(R = R.version.string, lavaan = as.character(packageVersion("lavaan")), psych = as.character(packageVersion("psych"))),
  n = nrow(DAT), seed = SEED, splits = N_SPLITS, candidates = length(SUBSETS), choices = CHOICES,
  statements = as.list(SRC[ST]), benchmarks = lapply(BENCH, lab),
  reference_check_full_sample = ref_chk, screen_check = SCHK,
  primary_split = list(search_n = length(SPLITS[[1]]$search), test_n = length(SPLITS[[1]]$test),
                       search_by_stratum = as.list(table(strata[SPLITS[[1]]$search])), test_by_stratum = as.list(table(strata[SPLITS[[1]]$test])),
                       walk = WALK[split == 1], best = TOPN, pareto_front = FRONT, bootstrap = BS),
  table = TAB[, ..SMALL],
  material_in_primary_test_half = TAB[sample == "test" & material_separation == TRUE & !grepl("benchmark", role), .(id, role, statements, version, r_obs, d_r_obs, omega)],
  stability = list(summary = OPT, subsets = STAB_SUB, statements = STAB_ST),
  failures = list(phase3_models_failed = VE[ok == FALSE, .(id, sample, version, block, error)],
                  phase3_models_improper = VE[ok == TRUE & nzchar(improper), unique(.SD), .SDcols = c("id", "sample", "version", "block", "improper")],
                  h1_not_met_or_improper = TAB[!(h1_rule_met %in% TRUE) | nzchar(trimws(improper_h1)), .(id, role, sample, h1_rule_met, improper_h1)],
                  search_walk_failures = WALK[, lapply(.SD, sum), by = .(version, objective), .SDcols = c("checked", "fail_h1_rule", "fail_omega", "fail_improper_or_orientation_h1", "fail_validity_model")]),
  run_seconds = list(search = round(secs_search), per_split = SPLIT_SECS, evaluation = round(secs_eval), bootstrap = round(secs_boot),
                     total = round(as.numeric(difftime(Sys.time(), T0, units = "secs")))), cores = CORES)
write_json(rnd(out), out21("21_exploratory_subsets.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
sf2b_log("21 done in ", out$run_seconds$total, " s; outputs ", out21("21_*"))
print(num(copy(TAB))[, ..SMALL])
print(OPT); print(STAB_SUB[times_chosen > 1]); print(BS); print(SCHK)
