# 14: which statements drive the 8's overlap with FFMQ Observing, NCS-6 and BIDR SDE.
# EXPLORATORY, NOT REGISTERED.
#
# Not part of the filed registration (osf.io/d4a5u) and not in its list of exploratory
# analyses. Specified on 15 November 2026, after the confirmatory real-mode run and
# analyses 11-13 had been seen, so every number here is hypothesis-generating only. No
# registered cell, test or output is changed.
#
# Why: the 8's latent r is .865 with FFMQ Observing, .453 with NCS-6 and .484 with BIDR
# Self-Deceptive Enhancement (04). 12 found three-quarters of the 8's variance shared with
# Observing at construct level, with no duplicated item pairs; 13 found a general
# self-description factor carrying most of the NCS-6 and BIDR SDE correlations. This
# asks which statements carry those overlaps.
#
# Statements: the 8, A9 and sfc_24 (the 10's two), and the three unscored reversed
# statements (N6, sf6_44, sfc_31): all 13. Named by id only.
#   (a) Total: each statement as a single ordinal indicator (a single-indicator factor
#       with no residual variance, 04's treatment of single items) beside each measure's
#       factor(s), the measure's own model as in 04 (FFMQ five facets, BIDR both scales,
#       BFI-2 three facets, Brief MAIA-2 eight subscales). Its latent r with every
#       battery measure; the three targets are the focus, the other 25 feed (c).
#   (b) Net of the 8's content factor: the H1 model (the 8 with its agreement factor, as
#       registered) plus the measure, and the statement's unique part as a factor of its
#       own (U, loading 1, variance 1, the statement's residual fixed at 0: under theta
#       this is the same model as residual variance 1), uncorrelated with S and the
#       agreement factor and correlated with the measure's factor(s). A statement outside
#       the 8 is added to the model as a ninth statement (free loading on S, loading 1 on
#       the agreement factor, as in 03's 13-statement model). Reported: r(U, measure),
#       the statement's specific association, and u * r(U, measure), the part of the
#       statement's latent-response correlation with the measure not carried by S
#       (u = U's standardized loading).
#   (c) The general factor: 13's Schmid-Leiman G defined without the 8 (13's
#       "sl_without_8", rebuilt from 13_latent_matrix.csv), each statement extended onto
#       it by Dwyer's (1937) method from its 28 latent correlations in (a). Cheap: (a)
#       already holds the correlations. Also the statement's correlation with each
#       target net of G.
#   (d) Leave-one-out: each of the 8 dropped in turn; the 7 refitted as H1 with its
#       agreement factor (fit, the largest residual and the local rule, omega, alpha) and
#       with each target's measure (04's model with the 7): its latent r. The 8 itself is
#       refitted in the same code as the reference.
#   Ranking: the 8 by how much the overlaps fall when each is dropped (r_8 - r_7, the
#   mean over the three targets); all 13 by their keyed specific association.
#
# Settings as registered (Section 5): ordinal items, WLSMV, theta, missing = "pairwise",
# conditional.x = FALSE, the 8's agreement factor as registered. Brief MAIA-2 enters with
# decision 46's bound on maia_12, as in 13. CHOICES below are written to the JSON.
#
# Writes output/<tag>/14_*.csv and 14_exploratory_items.json (aggregates only).
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/14_exploratory_items.R

# Set Up ---------
## Load libraries ---------
packages <- c("jsonlite", "data.table", "lavaan", "psych", "GPArotation")
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
# 10-13 do): A, D, M, S, DAT, score_items; nothing is fitted.
ex04 <- parse("study2b/analysis/04_validity.R", keep.source = FALSE)
first_fit <- which(vapply(ex04, function(e) grepl("^sf2b_log\\(\"fitting", deparse(e)[1]), TRUE))
stopifnot(length(first_fit) == 1)
for (i in seq_len(first_fit - 1)) eval(ex04[[i]], envir = globalenv())
stopifnot(exists("DAT"), exists("score_items"), !exists("RES"))
ex10 <- parse("study2b/analysis/10_deviation_maia.R", keep.source = FALSE)
is_def <- function(e, nm) is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(nm))
eval(Filter(function(e) is_def(e, "improper_of"), as.list(ex10))[[1]], envir = globalenv())
improper_list <- function(f) { x <- improper_of(f)$list; x[!grepl("^std loading (L_[[:alnum:]_]+|grp) =~ ", x)] }

I8 <- D$sf8
ST <- c(I8, "sf2b_08", "sf2b_12", "sf2b_05", "sf2b_11", "sf2b_13")
stopifnot(setequal(ST, D$items))
SRC <- D$source; REV <- D$reversed
key <- function(it) ifelse(it %in% REV, -1, 1)
role_of <- function(it) ifelse(it %in% I8, "the 8", ifelse(it %in% D$sf10, "the 10", "unscored reversed"))
TARGETS <- c("ffmq_observing", "ncs6", "bidr_sde")
BLOCKS <- list(ffmq = SF2B_MEASURES$ffmq, sbi = "sbi_moment", bfi = c("bfi_aesthetic", "bfi_intellectual", "bfi_creative"),
               eq = "eq_decentering", samd = "samd_md", ncs = "ncs6", bidr = c("bidr_sde", "bidr_im"), bsss = "bsss8",
               spsq = "spsq_sd", asi = "asi_pc", maia = SF2B_MEASURES$maia, pss = "pss4", phq = "phq4_total", lifesat = "lifesat")
TARGET_BLOCK <- c(ffmq_observing = "ffmq", ncs6 = "ncs", bidr_sde = "bidr")
MAIA_BOUND <- 9   # decision 46 under theta: l^2 <= 9, std loading <= .949
OMEGA_TOL <- .03  # the registered H1 tolerance for a change in omega, used as the yardstick for "keeping omega"

CHOICES <- list(
  statements = "All 13: the 8, A9 (sf2b_08), sfc_24 (sf2b_12) and the unscored reversed N6 (sf2b_05), sf6_44 (sf2b_11), sfc_31 (sf2b_13). Ids only; no wording.",
  keyed_sign = "Correlations and loadings are reported raw and keyed (the reversed statements' signs flipped), so a keyed positive value means the statement, scored as the scale scores it, goes with the measure.",
  total = "Each statement as a single-indicator ordinal factor (loading 1, residual 0, variance 1: 04's treatment of single items) beside the measure's factor(s) in the measure's own model as 04 builds it (FFMQ five facets; BIDR SDE and IM; BFI-2 three facets; Brief MAIA-2 eight subscales with decision 46's bound on maia_12). No agreement factor: one statement cannot identify it, so these correlations carry the statement's acquiescence.",
  net_of_content = "The H1 model (the 8 + agreement factor, registered) plus the measure's model (04's syntax, factors uncorrelated with the agreement factor), plus U, the statement's unique part (U =~ 1*statement, statement residual fixed at 0, U variance 1; the same model as theta's residual variance 1), U uncorrelated with S and the agreement factor and free with every factor of the measure. A statement outside the 8 enters as a ninth statement (free loading on S, loading 1 on the agreement factor, as 03's 13-statement model). Reported: r(U, target) with its CI (the specific association), u*r(U, target) (the semipartial part of the statement's latent-response correlation not carried by S), and r(S, target) in the same model.",
  general_factor = "13's Schmid-Leiman general factor defined without the 8 (sl_without_8: minres, oblimin, k from 13's parallel analysis), rebuilt from 13_latent_matrix.csv (4 decimals) with 13's smoothing; each statement extended by Dwyer's (1937) least-squares extension of the oblique first-order solution through the second-order loadings, exactly 13's sl_extend(), from its 28 latent correlations in (a). The rebuild is checked against 13_general_loadings.csv. Statement partial r net of G: (r - g_i g_F) / sqrt((1 - g_i^2)(1 - g_F^2)), g_F the target's loading on the same G. Point estimates only (no SEs): the matrix is assembled from separate models.",
  leave_one_out = "Each of the 8 dropped in turn; the 7 with its agreement factor (registered form) through sf2b_h1 (fit, largest residual, local rule, close-fit test, content omega) and alpha as 03 computes it; and the 7 in 04's model of each target's measure (FFMQ five facets, NCS-6, BIDR both scales) for its latent r. The full 8 refitted in the same code as the reference.",
  agreement_sensitivity = "The registered model keeps the agreement factor uncorrelated with every validity factor. 04 found it correlates .35 with FFMQ Observing when freed (04_agreement_factor_validity.csv), so acquiescence shared by the 8's unreversed statements and Observing's (all positively worded) items can only pass through S or the statements' unique parts. The net-of-content and leave-one-out models are therefore also fitted with the agreement factor free to correlate with the measure's factor(s) (agreement_free = TRUE); the registered form stays primary and the ranking uses it.",
  ranking = "The 8 ranked by the mean over the three targets of r_8 - r_7 (how much the overlap falls without the statement). A candidate for removal or rewriting: the mean fall is at least .01, the 7 still meets the registered H1 rule, and its omega is within .03 (the registered H1 omega tolerance) of the 8's. All 13 also ranked by the mean keyed semipartial specific association over the three targets.",
  improper = "10's improper_of() on every fitted model (std loading >= 1, negative residual or factor variances, factor correlations >= 1, non-positive-definite factor covariance); single-indicator factors' loadings of 1 are by construction and not listed. A model that fails, is not oriented as keyed, or is improper is reported and its estimates are flagged.")

## Syntax ---------
fac_of <- function(sc) if (sc == "lifesat") "L_lifesat" else sc
block_lines <- function(b, ri, ri_free = FALSE) {
  if (b == "lifesat") return(list(syn = c("L_lifesat =~ 1*lifesat", "lifesat ~~ 0*lifesat", "L_lifesat ~~ 1*L_lifesat", if (ri) "RI ~~ 0*L_lifesat"),
                                  ord = "lifesat", scores = character()))
  sc <- BLOCKS[[b]]
  x <- sf2b_validity_syntax(sc, D, ri = ri, ri_free = ri_free)
  x <- x[!grepl("^S ~~", x)]
  if ("maia_notworrying" %in% sc) x <- c(x, "maia_notworrying =~ l_maia12*maia_12", sprintf("l_maia12^2 < %s", MAIA_BOUND))
  list(syn = x, ord = unlist(score_items[sc]), scores = sc)
}
orient_validity <- function(f, scores) {
  if (!length(scores)) return(character())
  L <- lavInspect(f, "std")$lambda; its <- sf2b_score_items(D); rv <- unlist(lapply(D$validity, `[[`, "reverse"))
  Filter(function(sc) sum(L[its[[sc]], sc] * ifelse(its[[sc]] %in% rv, -1, 1)) <= 0, scores)
}
model_total <- function(st, b) {
  bl <- block_lines(b, ri = FALSE); fs <- vapply(BLOCKS[[b]], fac_of, "")
  list(syn = c(paste0("L_", st, " =~ 1*", st), paste0(st, " ~~ 0*", st), paste0("L_", st, " ~~ 1*L_", st), bl$syn,
               paste0("L_", st, " ~~ r_", BLOCKS[[b]], "*", fs)),
       ord = unique(c(st, bl$ord)), scores = bl$scores)
}
model_net <- function(st, b, ri_free = FALSE) {
  items <- if (st %in% I8) I8 else c(I8, st)
  bl <- block_lines(b, ri = TRUE, ri_free = ri_free); sc <- BLOCKS[[b]]
  list(syn = c(sf2b_scale_syntax(items), bl$syn, paste0("S ~~ r_", sc, "*", sc),
               paste0("U_", st, " =~ 1*", st), paste0(st, " ~~ 0*", st), paste0("U_", st, " ~~ 1*U_", st),
               paste0("U_", st, " ~~ 0*S"), paste0("U_", st, " ~~ 0*RI"), paste0("U_", st, " ~~ e_", sc, "*", sc)),
       ord = unique(c(items, bl$ord)), items = items, scores = sc)
}
# S takes its sign from its first statement (lavaan's NA* default): an unreversed statement
# is put first, so dropping N4 (the first) does not leave S oriented on reversed sf6_36.
model_loo <- function(items, b, ri_free = FALSE) {
  items <- items[order(items %in% REV)]
  bl <- block_lines(b, ri = TRUE, ri_free = ri_free); sc <- BLOCKS[[b]]
  list(syn = c(sf2b_scale_syntax(items), bl$syn, paste0("S ~~ r_", sc, "*", sc)), ord = unique(c(items, bl$ord)), items = items, scores = sc)
}

## Tasks ---------
TASKS <- list()
add <- function(...) TASKS[[length(TASKS) + 1L]] <<- list(...)
for (st in ST) for (b in names(BLOCKS)) add(kind = "total", st = st, block = b)
for (rf in c(FALSE, TRUE)) for (st in ST) for (tg in TARGETS) add(kind = "net", st = st, block = TARGET_BLOCK[[tg]], ri_free = rf)
for (dr in c(NA, I8)) {
  add(kind = "h1", drop = dr)
  for (rf in c(FALSE, TRUE)) for (tg in TARGETS) add(kind = "loo", drop = dr, block = TARGET_BLOCK[[tg]], ri_free = rf)
}
run_task <- function(t) {
  t0 <- Sys.time()
  base <- list(task = t, ok = FALSE, error = "", improper = "")
  if (t$kind == "h1") {
    items <- setdiff(I8, t$drop); items <- items[order(items %in% REV)]
    h <- sf2b_try(sf2b_h1(X, items, REV))
    if (inherits(h, "error") || !isTRUE(h$converged)) { base$error <- if (inherits(h, "error")) conditionMessage(h) else h$error; return(base) }
    return(c(base[c("task")], list(ok = TRUE, error = "", improper = paste(improper_list(h$fit), collapse = "; "),
             stats = h$stats, max_residual = h$max_residual, max_residual_pair = h$max_residual_pair, d_omega = h$d_omega, d_lambda = h$d_lambda,
             local_ok = h$local_ok, cfi_ok = h$cfi_ok, close_ok = h$close_ok, closefit_p = h$closefit_p, supported = h$supported,
             omega = h$omega, loadings = h$loadings, ri_share = h$ri_share, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))))
  }
  m <- switch(t$kind, total = model_total(t$st, t$block), net = model_net(t$st, t$block, t$ri_free), loo = model_loo(setdiff(I8, t$drop), t$block, t$ri_free))
  f <- sf2b_try(sf2b_sem(m$syn, DAT, ordered = m$ord))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (sf2b_failed(f)) { base$error <- sf2b_msg(f); base$secs <- secs; return(base) }
  bad <- if (t$kind == "total") orient_validity(f, m$scores) else sf2b_orientation(f, m$items, m$scores, D)
  pe <- parameterEstimates(f); ss <- standardizedSolution(f)
  est <- pe[grepl("^(r|e)_", pe$label), c("label", "est", "se", "ci.lower", "ci.upper", "pvalue")]
  if (isTRUE(t$ri_free)) { x <- ss[grepl("^ri_", ss$label), c("label", "est.std", "se", "ci.lower", "ci.upper", "pvalue")]; names(x)[2] <- "est"; est <- rbind(est, x) }
  if (t$kind == "net") {
    e <- ss[grepl("^e_", ss$label), c("label", "est.std", "se", "ci.lower", "ci.upper", "pvalue")]; names(e)[2] <- "est"
    est <- rbind(est[grepl("^(r|ri)_", est$label), ], e)
    std <- lavInspect(f, "std")$lambda
    extra <- list(u = std[t$st, paste0("U_", t$st)], lam = std[t$st, "S"], agr = std[t$st, "RI"])
  } else extra <- list()
  c(list(task = t, ok = !length(bad), error = if (length(bad)) paste("not oriented as keyed:", paste(bad, collapse = ", ")) else "",
         improper = paste(improper_list(f), collapse = "; "), est = est, fit = sf2b_fit_stats(f), secs = secs), extra)
}
sf2b_log("14: ", length(TASKS), " models on ", SF2B$cores, " workers")
t1 <- Sys.time()
RES <- sf2b_par_lapply(TASKS, run_task, data = list(DAT = DAT, X = A$scale_raw, D = D, score_items = score_items, I8 = I8, REV = REV,
                                                     BLOCKS = BLOCKS, MAIA_BOUND = MAIA_BOUND, fac_of = fac_of, block_lines = block_lines,
                                                     orient_validity = orient_validity, model_total = model_total, model_net = model_net,
                                                     model_loo = model_loo, improper_of = improper_of, improper_list = improper_list))
secs_models <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
tid <- function(t) paste(t$kind, t$st %||% (if (is.null(t$drop) || is.na(t$drop)) "none dropped" else paste("drop", t$drop)), t$block %||% "",
                         if (isTRUE(t$ri_free)) "agreement free" else "", sep = " | ")
for (r in Filter(function(r) !r$ok, RES)) sf2b_log("FAILED: ", tid(r$task), ": ", r$error)
for (r in Filter(function(r) r$ok && nzchar(r$improper), RES)) sf2b_log("IMPROPER: ", tid(r$task), ": ", r$improper)
fits <- rbindlist(lapply(RES, function(r) c(list(kind = r$task$kind, statement = r$task$st %||% NA_character_,
                                                  dropped = if (is.null(r$task$drop)) NA_character_ else r$task$drop,
                                                  block = r$task$block %||% NA_character_, agreement_free = isTRUE(r$task$ri_free), ok = r$ok, error = r$error, improper = r$improper %||% "",
                                                  seconds = round(r$secs %||% NA_real_, 1)),
                                             if (r$ok && !is.null(r$fit)) r$fit else if (r$ok) r$stats)), fill = TRUE)
fwrite(fits, sf2b_out("14_model_fit.csv"))
num <- function(DT, k = 3) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, function(x) round(x, k)), .SDcols = cols]; DT }

## (a) Total ---------
TOT <- rbindlist(lapply(Filter(function(r) r$ok && r$task$kind == "total", RES), function(r)
  data.table(statement = r$task$st, block = r$task$block, as.data.table(r$est), improper = r$improper)))
TOT[, measure := sub("^r_", "", label)]
TOT[, `:=`(source = SRC[statement], role = role_of(statement), k = key(statement))]
TOT[, `:=`(r_keyed = est * k, lower_keyed = fifelse(k > 0, ci.lower, -ci.upper), upper_keyed = fifelse(k > 0, ci.upper, -ci.lower))]
TOT_out <- TOT[, .(statement, source, role, measure, block, r_raw = est, se, r_keyed, lower_keyed, upper_keyed, p = signif(pvalue, 3), improper)]
TOT_out <- TOT_out[order(match(statement, ST), match(measure, unlist(BLOCKS)))]
fwrite(num(copy(TOT_out)), sf2b_out("14_item_correlations.csv"))

## (b) Net of the 8's content factor ---------
NET <- rbindlist(lapply(Filter(function(r) r$ok && r$task$kind == "net", RES), function(r) {
  tg <- names(TARGET_BLOCK)[TARGET_BLOCK == r$task$block]
  e <- r$est[r$est$label %in% paste0("e_", tg), ]; s <- r$est[r$est$label %in% paste0("r_", tg), ]
  a <- r$est[r$est$label %in% paste0("ri_", tg), ]
  data.table(statement = r$task$st, target = tg, agreement_free = isTRUE(r$task$ri_free), r_agreement_target = if (nrow(a)) a$est else NA_real_,
             r_S_target = s$est, specific_r = e$est, specific_lower = e$ci.lower, specific_upper = e$ci.upper,
             specific_p = e$pvalue, u = r$u, content_loading = r$lam, agreement_loading = r$agr, improper = r$improper)
}))
NET[, k := key(statement)]
NET[, `:=`(specific_r_keyed = specific_r * k, specific_lower_keyed = fifelse(k > 0, specific_lower, -specific_upper),
           specific_upper_keyed = fifelse(k > 0, specific_upper, -specific_lower), semipartial_keyed = u * specific_r * k,
           via_S_keyed = content_loading * k * r_S_target, content_loading_keyed = content_loading * k)]
NET[, model_total_keyed := via_S_keyed + semipartial_keyed]
NET <- merge(NET, TOT[measure %in% TARGETS, .(statement, target = measure, single_indicator_r_keyed = r_keyed)], by = c("statement", "target"), all.x = TRUE)
NET[, `:=`(source = SRC[statement], role = role_of(statement))]
NET_out <- NET[order(agreement_free, match(target, TARGETS), match(statement, ST)),
               .(agreement_free, target, statement, source, role, single_indicator_r_keyed, r_agreement_target, content_loading_keyed, agreement_loading, r_S_target,
                 via_S_keyed, specific_r_keyed, specific_lower_keyed, specific_upper_keyed, specific_p = signif(specific_p, 3),
                 semipartial_keyed, model_total_keyed, improper)]
fwrite(num(copy(NET_out)), sf2b_out("14_item_specific.csv"))

## (c) The general factor (13's G without the 8), extended to the statements ---------
LM <- fread(sf2b_out("13_latent_matrix.csv")); RL <- as.matrix(LM[, -1, with = FALSE]); rownames(RL) <- LM$var
J13 <- fromJSON(sf2b_out("13_exploratory_general_factor.json"))
BAT <- J13$matrix$battery; K <- as.integer(J13$parallel_analysis$group_factors_used); N <- nrow(DAT)
Rm <- RL[BAT, BAT]; Rm <- (Rm + t(Rm)) / 2
smoothed <- min(eigen(Rm, symmetric = TRUE, only.values = TRUE)$values) <= 0
if (smoothed) Rm <- suppressWarnings(cor.smooth(Rm))
sl_fit <- function(R, k) {   # 13's
  s <- schmid(R, nfactors = k, fm = "minres", rotate = "oblimin", n.obs = N)
  g <- setNames(s$sl[, "g"], rownames(s$sl))
  flip <- if (isTRUE(g["sf8"] < 0) || (!"sf8" %in% names(g) && sum(g) < 0)) -1 else 1
  list(s = s, g = g * flip, flip = flip, k = k)
}
sl_extend <- function(sl, Roe) {   # 13's
  P <- unclass(sl$s$oblique); Phi <- sl$s$phi; Sg <- P %*% Phi
  pe <- Roe[, rownames(P), drop = FALSE] %*% Sg %*% solve(crossprod(Sg))
  gam <- as.numeric(unclass(sl$s$gloading))
  drop(pe %*% gam) * sl$flip
}
SL <- suppressWarnings(sl_fit(Rm, K))
B28 <- setdiff(BAT, "sf8")
SL28 <- suppressWarnings(sl_fit(Rm[B28, B28], K))
if (cor(SL28$g, SL$g[B28]) < 0) { SL28$g <- -SL28$g; SL28$flip <- -SL28$flip }
g8_ext <- unname(sl_extend(SL28, Rm["sf8", B28, drop = FALSE]))[1]
L13 <- fread(sf2b_out("13_general_loadings.csv"))
rebuild_check <- list(max_abs_diff_sl_without_8 = max(abs(c(SL28$g, sf8 = g8_ext)[L13$measure] - L13$sl_g_without_8)),
                      max_abs_diff_sl = max(abs(SL$g[L13$measure] - L13$sl_g)), g8_extended = g8_ext, smoothed = smoothed, k = K)
sf2b_log("13's G rebuilt: max |difference| ", signif(rebuild_check$max_abs_diff_sl_without_8, 2), " (sl_without_8)")
stopifnot(rebuild_check$max_abs_diff_sl_without_8 < .01)
RO <- dcast(TOT[measure %in% B28], statement ~ measure, value.var = "est")
for (b in setdiff(B28, names(RO))) RO[, (b) := NA_real_]
complete28 <- RO[, statement][complete.cases(RO[, ..B28])]
Roe <- as.matrix(RO[statement %in% complete28, ..B28]); rownames(Roe) <- RO[statement %in% complete28, statement]
g_raw <- sl_extend(SL28, Roe)
GEN <- data.table(statement = names(g_raw), source = SRC[names(g_raw)], role = role_of(names(g_raw)), g_loading_raw = g_raw,
                  g_loading_keyed = g_raw * key(names(g_raw)))
H1REF <- fread(sf2b_out("03_loadings.csv"))[form == "13", .(statement = item, content_loading_13 = keyed_content, agreement_loading_13 = agreement)]
GEN <- merge(GEN, H1REF, by = "statement", all.x = TRUE)
for (tg in TARGETS) {
  gF <- SL28$g[[tg]]
  x <- TOT[measure == tg, .(statement, r_keyed)]
  GEN <- merge(GEN, x[, .(statement, rk = r_keyed)], by = "statement", all.x = TRUE)
  GEN[, (paste0("r_", tg)) := rk]
  GEN[, (paste0("net_of_G_", tg)) := (rk - g_loading_keyed * gF) / sqrt((1 - g_loading_keyed^2) * (1 - gF^2))]
  GEN[, (paste0("share_via_G_", tg)) := (g_loading_keyed * gF) / rk]
  GEN[, rk := NULL]
}
GEN <- GEN[order(match(statement, ST))]
fwrite(num(copy(GEN)), sf2b_out("14_item_general.csv"))
target_g <- SL28$g[TARGETS]

## (d) Leave-one-out ---------
dkey <- function(t) if (is.null(t$drop) || is.na(t$drop)) "none" else t$drop
H1T <- rbindlist(lapply(Filter(function(r) r$task$kind == "h1", RES), function(r) {
  if (!r$ok) return(data.table(dropped = dkey(r$task), ok = FALSE, error = r$error))
  data.table(dropped = dkey(r$task), ok = TRUE, error = "", n = r$stats$ntotal, cfi = r$stats$cfi.scaled, rmsea = r$stats$rmsea.scaled,
             rmsea_lower = r$stats$rmsea.ci.lower.scaled, rmsea_upper = r$stats$rmsea.ci.upper.scaled, srmr = r$stats$srmr,
             max_residual = r$max_residual, max_residual_pair = r$max_residual_pair, d_omega_if_freed = r$d_omega, d_lambda_if_freed = r$d_lambda,
             cfi_ok = r$cfi_ok, local_ok = r$local_ok, close_ok = r$close_ok, h1_rule_met = r$supported, omega = r$omega,
             agreement_share = r$ri_share, improper = r$improper)
}), fill = TRUE)
alpha_of <- function(items) { Z <- A$scale_raw[, items]; rv <- intersect(items, REV); Z[, rv] <- 7 - Z[, rv]   # 03's
  suppressWarnings(psych::alpha(Z[complete.cases(Z), ], check.keys = FALSE, warnings = FALSE)$total$raw_alpha) }
H1T[, alpha := vapply(dropped, function(d) alpha_of(setdiff(I8, d)), 0)]
LOOR <- rbindlist(lapply(Filter(function(r) r$task$kind == "loo", RES), function(r) {
  if (!r$ok) return(data.table(dropped = dkey(r$task), block = r$task$block, agreement_free = r$task$ri_free, ok = FALSE, error = r$error))
  x <- r$est[r$est$label %in% paste0("r_", TARGETS), ]
  data.table(dropped = dkey(r$task), block = r$task$block, agreement_free = r$task$ri_free, ok = TRUE, error = "", target = sub("^r_", "", x$label),
             r = x$est, lower = x$ci.lower, upper = x$ci.upper, improper = r$improper)
}), fill = TRUE)
LW <- dcast(LOOR[ok == TRUE & agreement_free == FALSE], dropped ~ target, value.var = c("r", "lower", "upper"))
LWf <- dcast(LOOR[ok == TRUE & agreement_free == TRUE], dropped ~ target, value.var = "r")
setnames(LWf, TARGETS, paste0("r_", TARGETS, "_agreement_free"))
LOO <- merge(merge(H1T, LW, by = "dropped", all.x = TRUE, sort = FALSE), LWf, by = "dropped", all.x = TRUE, sort = FALSE)
ref <- LOO[dropped == "none"]
for (tg in TARGETS) LOO[, (paste0("fall_", tg)) := ref[[paste0("r_", tg)]] - get(paste0("r_", tg))]
for (tg in TARGETS) LOO[, (paste0("fall_", tg, "_agreement_free")) := ref[[paste0("r_", tg, "_agreement_free")]] - get(paste0("r_", tg, "_agreement_free"))]
LOO[, mean_fall := rowMeans(.SD), .SDcols = paste0("fall_", TARGETS)]
LOO[, mean_fall_agreement_free := rowMeans(.SD), .SDcols = paste0("fall_", TARGETS, "_agreement_free")]
LOO[, `:=`(omega_change = omega - ref$omega, alpha_change = alpha - ref$alpha)]
LOO[, source := fifelse(dropped == "none", "(the 8)", SRC[dropped])]
setcolorder(LOO, c("dropped", "source"))
LOO <- LOO[order(dropped != "none", -mean_fall)]
fwrite(num(copy(LOO)), sf2b_out("14_leave_one_out.csv"))
chk04 <- fread(sf2b_out("04_latent_validity.csv"))[form == "sf8" & sample == "all" & ri_free == FALSE & score %in% TARGETS, .(score, r_04 = est)]
chk04[, r_here := unlist(ref[, paste0("r_", score), with = FALSE])]

## Ranking ---------
spec_w <- dcast(NET[agreement_free == FALSE], statement ~ target, value.var = "semipartial_keyed")
setnames(spec_w, TARGETS, paste0("semipartial_", TARGETS))
spec_w[, mean_semipartial := rowMeans(.SD), .SDcols = paste0("semipartial_", TARGETS)]
er <- dcast(NET[agreement_free == FALSE], statement ~ target, value.var = "specific_r_keyed"); setnames(er, TARGETS, paste0("specific_r_", TARGETS))
spf <- dcast(NET[agreement_free == TRUE], statement ~ target, value.var = "semipartial_keyed"); setnames(spf, TARGETS, paste0("semipartial_", TARGETS, "_agreement_free"))
spf[, mean_semipartial_agreement_free := rowMeans(.SD), .SDcols = paste0("semipartial_", TARGETS, "_agreement_free")]
er <- merge(er, spf, by = "statement")
RK <- merge(merge(spec_w, er, by = "statement"), GEN[, .(statement, g_loading_keyed, content_loading_13)], by = "statement", all.x = TRUE)
RK <- merge(RK, LOO[dropped != "none", c("dropped", paste0("fall_", TARGETS), "mean_fall", "mean_fall_agreement_free", "omega", "omega_change", "h1_rule_met"), with = FALSE],
            by.x = "statement", by.y = "dropped", all.x = TRUE)
RK[, `:=`(source = SRC[statement], role = role_of(statement))]
RK[, rank_specific_all13 := frank(-mean_semipartial, ties.method = "min")]
RK[statement %in% I8, rank_loo_fall := frank(-mean_fall, ties.method = "min")]
RK[, candidate := statement %in% I8 & mean_fall >= .01 & h1_rule_met %in% TRUE & omega_change >= -OMEGA_TOL]
setcolorder(RK, c("statement", "source", "role", "rank_loo_fall", "mean_fall", paste0("fall_", TARGETS), "omega", "omega_change", "h1_rule_met",
                  "candidate", "rank_specific_all13", "mean_semipartial"))
RK <- RK[order(is.na(rank_loo_fall), rank_loo_fall, rank_specific_all13)]
fwrite(num(copy(RK)), sf2b_out("14_item_ranking.csv"))

## Report ---------
rnd <- function(x) if (is.data.frame(x)) x else if (is.list(x)) lapply(x, rnd) else if (is.numeric(x)) signif(x, 4) else x
out <- list(
  analysis = "14: which statements drive the 8's overlap with FFMQ Observing, NCS-6 and BIDR SDE. EXPLORATORY, NOT REGISTERED (osf.io/d4a5u does not include it); specified after the confirmatory run and analyses 11-13, with their results seen.",
  mode = SF2B$mode, tag = SF2B$tag, written = format(Sys.time(), "%Y-%m-%d %H:%M"),
  software = list(R = R.version.string, lavaan = as.character(packageVersion("lavaan")), psych = as.character(packageVersion("psych")),
                  GPArotation = as.character(packageVersion("GPArotation"))),
  n = nrow(DAT), choices = CHOICES, statements = as.list(SRC[ST]),
  models = list(fitted = length(RES), failed = fits[ok == FALSE, .(kind, statement, dropped, block, error)],
                improper = fits[ok == TRUE & nzchar(improper), .(kind, statement, dropped, block, improper)]),
  reference_8_against_04 = chk04,
  total_targets = TOT_out[measure %in% TARGETS], net_of_content = NET_out,
  general_factor = list(rebuild_check = rebuild_check, target_g_loadings = as.list(target_g), statements = GEN,
                        statements_without_all_28_cells = setdiff(ST, complete28)),
  leave_one_out = LOO, ranking = RK, candidates = RK[candidate == TRUE, .(statement, source, mean_fall, omega, omega_change)],
  run_seconds = list(models = round(secs_models), total = round(as.numeric(difftime(Sys.time(), T0, units = "secs")))), cores = SF2B$cores)
write_json(rnd(out), sf2b_out("14_exploratory_items.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
sf2b_log("14 done in ", out$run_seconds$total, " s; outputs ", sf2b_out("14_*"))
print(chk04)
print(num(copy(NET_out))[, .(agreement_free, target, source, single_indicator_r_keyed, content_loading_keyed, specific_r_keyed, specific_lower_keyed, specific_upper_keyed, semipartial_keyed)])
print(num(copy(GEN))[, c("source", "g_loading_keyed", "content_loading_13", paste0("net_of_G_", TARGETS)), with = FALSE])
print(num(copy(LOO))[, .(source, cfi, rmsea, srmr, max_residual, h1_rule_met, omega, alpha, r_ffmq_observing, r_ncs6, r_bidr_sde, mean_fall, r_ffmq_observing_agreement_free, r_ncs6_agreement_free, r_bidr_sde_agreement_free, mean_fall_agreement_free)])
print(num(copy(RK))[, .(source, rank_loo_fall, mean_fall, mean_fall_agreement_free, omega_change, candidate, rank_specific_all13, mean_semipartial, mean_semipartial_agreement_free, g_loading_keyed)])
