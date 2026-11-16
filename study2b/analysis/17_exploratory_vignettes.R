# 17: the vignettes. EXPLORATORY, NOT REGISTERED.
#
# Not part of the filed registration (osf.io/d4a5u). Section 5 lists the vignettes among
# its exploratory analyses only as H8 (invariance across telling the modes apart, run in
# 08) and 08's descriptives; everything here goes beyond that and was specified on 15
# November 2026, after 08's results had been seen (455 of 720 tell the modes apart;
# doing vignettes 2 and 7 and sensing vignette 4 near the midpoint; the self-rating r .28
# with the 8; H8 invariance held through loadings). Every number is hypothesis-generating
# only, and no registered cell, test or output is changed.
#
#   1. Item analysis of the nine: mean, SD and the distribution of answers (1 = Doing
#      mode ... 6 = Sensing mode); the share on the intended side (sensing >= 4, doing
#      <= 3); each vignette's correlation with the separation score recomputed without
#      it (raw, and keyed so that positive means it separates as intended); its
#      correlation with the observed 8 and with the self-rating. Vignettes 2, 4 and 7
#      flagged as ambiguous (named in advance from 08's means). Also alpha of the nine
#      keyed vignettes, with and without the three.
#   2. Re-scored without 2, 4 and 7 (sensing 3, 6, 8; doing 1, 5, 9): separation, the
#      registered cut (>= 1) and the share who tell the modes apart, against the
#      registered classification.
#   3. Moderation by telling the modes apart (yes/no, the registered definition, 455 vs
#      265): the 8's latent r with FFMQ Observing, SBI, BFI-2 Aesthetic Sensitivity,
#      Brief MAIA-2 Noticing, NCS-6, BIDR SDE, past-week pause days and the self-rating,
#      and the H6 difference (MIMIC), in two-group models with the parameter free vs
#      equal: the difference with its Wald CI, the scaled chi-square difference test,
#      and the smallest difference detectable at 80% power.
#   4. The self-rating against the 8 within each group.
#
# Settings as registered (Section 5): ordinal items, WLSMV, theta, missing = "pairwise",
# conditional.x = FALSE, the 8's agreement factor as registered (in each group: loadings
# 1, variance free, uncorrelated with every other factor). CHOICES below are in the JSON.
#
# Writes output/<tag>/17_*.csv and 17_exploratory_vignettes.json (aggregates only).
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/17_exploratory_vignettes.R

# Set Up ---------
## Load libraries ---------
packages <- c("jsonlite", "data.table", "lavaan", "psych")
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
ex04 <- parse("study2b/analysis/04_validity.R", keep.source = FALSE)
first_fit <- which(vapply(ex04, function(e) grepl("^sf2b_log\\(\"fitting", deparse(e)[1]), TRUE))
stopifnot(length(first_fit) == 1)
for (i in seq_len(first_fit - 1)) eval(ex04[[i]], envir = globalenv())
stopifnot(exists("DAT"), exists("score_items"), !exists("RES"))
ex10 <- parse("study2b/analysis/10_deviation_maia.R", keep.source = FALSE)
is_def <- function(e, nm) is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(nm))
eval(Filter(function(e) is_def(e, "improper_of"), as.list(ex10))[[1]], envir = globalenv())
not_improper <- function(x) grepl("^(group [[:alnum:]]+: )?std loading (L_[[:alnum:]_]+|grp) =~ ", x)
improper_list <- function(f) { x <- improper_of(f)$list; x[!not_improper(x)] }
# 10's improper_of() reads single-group matrices; the same checks, group by group.
improper_mg <- function(f) {
  if (lavInspect(f, "ngroups") == 1) return(improper_list(f))
  std <- lavInspect(f, "std"); est <- lavInspect(f, "est"); C <- lavInspect(f, "cor.lv"); CV <- lavInspect(f, "cov.lv")
  gl <- lavInspect(f, "group.label"); out <- character()
  for (g in seq_along(gl)) {
    p <- paste0("group ", gl[g], ": ")
    L <- std[[g]]$lambda; ix <- which(abs(L) >= 1, arr.ind = TRUE)
    if (nrow(ix)) out <- c(out, sprintf("%sstd loading %s =~ %s = %.3f", p, colnames(L)[ix[, 2]], rownames(L)[ix[, 1]], L[ix]))
    th <- diag(est[[g]]$theta); out <- c(out, if (any(th < 0)) sprintf("%snegative residual variance %s = %.3f", p, names(th)[th < 0], th[th < 0]))
    sth <- diag(std[[g]]$theta); out <- c(out, if (any(sth < 0)) sprintf("%snegative std residual variance %s = %.3f", p, names(sth)[sth < 0], sth[sth < 0]))
    ps <- diag(est[[g]]$psi); out <- c(out, if (any(ps < 0)) sprintf("%snegative factor variance %s = %.3f", p, names(ps)[ps < 0], ps[ps < 0]))
    Cg <- C[[g]]; C0 <- Cg; C0[lower.tri(C0, TRUE)] <- 0; cx <- which(abs(C0) >= 1, arr.ind = TRUE)
    if (nrow(cx)) out <- c(out, sprintf("%sfactor correlation %s ~~ %s = %.3f", p, rownames(Cg)[cx[, 1]], colnames(Cg)[cx[, 2]], C0[cx]))
    ev <- min(eigen(CV[[g]], symmetric = TRUE, only.values = TRUE)$values)
    if (ev <= 0) out <- c(out, sprintf("%sfactor covariance matrix not positive definite (smallest eigenvalue %.4f)", p, ev))
  }
  out[!not_improper(out)]
}
stopifnot(identical(M$participant, S$participant), nrow(DAT) == nrow(M))

I8 <- D$sf8; VK <- D$vkind; VIGS <- names(VK)
SENS <- VIGS[VK == "sensing"]; DOING <- VIGS[VK == "doing"]
AMBIG <- c("vig_2", "vig_4", "vig_7")
stopifnot(identical(SENS, c("vig_3", "vig_4", "vig_6", "vig_8")), identical(DOING, c("vig_1", "vig_2", "vig_5", "vig_7", "vig_9")),
          all(AMBIG %in% VIGS))
V <- A$vignettes
TA <- M$tells_apart
stopifnot(sum(TA == "yes", na.rm = TRUE) == 455, sum(TA == "no", na.rm = TRUE) == 265)
TARGETS <- c("ffmq_observing", "sbi_moment", "bfi_aesthetic", "maia_noticing", "ncs6", "bidr_sde", "crit_pause_days", "vig_self")
SINGLE <- c("crit_pause_days", "vig_self")
DAT$vig_self <- V[, "vig_self"]
DAT$grp <- fifelse(M$practice_group == "current", 1L, fifelse(M$practice_group == "never", 0L, NA_integer_))
DAT$.ta <- TA

CHOICES <- list(
  scale = "Each vignette is answered 1 (Doing mode) to 6 (Sensing mode). Separation (Section 4) = mean of the sensing vignettes - mean of the doing vignettes, at least three of each answered; tells the modes apart = separation >= 1.",
  ambiguous = "Vignettes 2 and 7 (doing) and 4 (sensing) flagged ambiguous in advance, from 08's means (3.82, 3.66, 3.66: near the scale's midpoint, 3.5).",
  rest_correlation = "Each vignette against the separation recomputed without it (the registered at-least-three rule on what remains). Keyed: doing vignettes sign-reversed, so a positive keyed correlation means the vignette separates the modes as the score does.",
  alpha = "Cronbach's alpha of the keyed vignettes (sensing as answered, doing 7 - x), complete cases; with all nine and without 2, 4 and 7. Descriptive: the vignettes are not built as a scale.",
  rescore = "Without 2, 4 and 7: sensing 3, 6, 8 and doing 1, 5, 9; the registered rule (at least three of each) then needs all six; tells apart = separation >= 1, the registered cut. Agreement with the registered classification: the cross-tab and Cohen's kappa.",
  moderation_groups = "Telling the modes apart as registered (yes 455, no 265). Two-group models, groups ordered no, yes.",
  moderation_model = "In each group: the 8's H1 model (S variance 1, loadings free, the agreement factor with loadings 1 and its variance free) and the target as one factor of its own items (04's syntax: variance 1, loadings started by keying, uncorrelated with the agreement factor), or as a single-indicator factor for pause days and the self-rating (04's treatment). Configural: loadings and thresholds free per group, latent-response residual variances fixed (theta), factor variances 1 in both groups, so S ~~ target is each group's latent correlation. Free: one correlation per group; equal: the same label in both. Each target alone (FFMQ Observing without the other facets, BIDR SDE without IM, Brief MAIA-2 Noticing without the other subscales, so maia_12 does not enter), to keep two-group models estimable; the pooled one-group fit of the same model is reported beside 04's cell.",
  sparse_categories = "Categories with fewer than 5 answers in either group merged toward the centre for that variable, in both groups (sf2b_merge_sparse, as H8), separately for each model's variables; merges are written out.",
  tests = "The difference r_yes - r_no from the free model (delta-method Wald CI and p), and lavTestLRT(equal, free): the scaled chi-square difference test for WLSMV.",
  detectable = "Smallest difference detectable at 80% power, two-sided .05: 2.80 x the SE of the difference in the free model (normal approximation). For observed correlations at n 455 and 265, a Fisher-z difference of 2.80 x sqrt(1/452 + 1/262) = .217, shown as the pair of correlations that differ by it around the pooled value.",
  h6 = "Current (1) vs never (0) practitioners as 07 codes them: two-group MIMIC (S ~ group, S's residual variance 1 in each group, so the coefficient is each group's latent d; agreement factor ~ 0*group), free vs equal; observed: 07's sf2b_kg within each group and OLS of the 8 on practice x tells apart.",
  observed = "Observed-score analogues: Pearson r of the observed 8 with the target's observed score within each group, and the Fisher-z test of their difference.",
  improper = "10's improper_of() checks, group by group for two-group models; single-indicator and covariate-phantom loadings of 1 are by construction and not listed.")

## Helpers ---------
num <- function(DT, k = 3) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, function(x) round(x, k)), .SDcols = cols]; DT }
rnd <- function(x) if (is.data.frame(x)) x else if (is.list(x)) lapply(x, rnd) else if (is.numeric(x)) signif(x, 4) else x
wilson <- function(k, n) prop.test(k, n, correct = FALSE)$conf.int
separation <- function(sens, doing, min_each = SF2B$min_vignettes_each) {
  s <- rowMeans(V[, sens, drop = FALSE], na.rm = TRUE) - rowMeans(V[, doing, drop = FALSE], na.rm = TRUE)
  s[rowSums(!is.na(V[, sens, drop = FALSE])) < min(min_each, length(sens)) | rowSums(!is.na(V[, doing, drop = FALSE])) < min(min_each, length(doing))] <- NA
  s
}
sep_all <- separation(SENS, DOING)
stopifnot(isTRUE(all.equal(unname(sep_all), unname(M$separation))))
r_ci <- function(x, y) sf2b_r_ci(x, y)

## 1. Item analysis ---------
sf2b_log("17: item analysis")
ITEM <- rbindlist(lapply(c(VIGS, "vig_self"), function(v) {
  x <- V[, v]; kind <- if (v == "vig_self") "self-rating" else VK[[v]]
  dist <- as.list(setNames(round(100 * tabulate(x[!is.na(x)], 6) / sum(!is.na(x)), 1), paste0("pct_", 1:6)))
  rest <- if (v == "vig_self") rep(NA_real_, length(x)) else separation(setdiff(SENS, v), setdiff(DOING, v))
  k <- if (identical(kind, "doing")) -1 else 1
  rr <- r_ci(x, rest); r8 <- r_ci(x, S$sf8); rs <- if (v == "vig_self") list(r = NA, lower = NA, upper = NA) else r_ci(x, V[, "vig_self"])
  c(list(vignette = v, kind = kind, flagged_ambiguous = v %in% AMBIG, n = sum(!is.na(x)), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE),
         distance_from_midpoint = abs(mean(x, na.rm = TRUE) - 3.5),
         pct_intended_side = if (kind == "sensing") 100 * mean(x >= 4, na.rm = TRUE) else if (kind == "doing") 100 * mean(x <= 3, na.rm = TRUE) else NA_real_),
    dist,
    list(r_rest_separation = rr$r, r_rest_lower = rr$lower, r_rest_upper = rr$upper, r_rest_keyed = k * rr$r,
         r_with_8 = r8$r, r_with_8_lower = r8$lower, r_with_8_upper = r8$upper,
         r_with_self = rs$r, r_with_self_lower = rs$lower, r_with_self_upper = rs$upper))
}), fill = TRUE)
keyed <- V[, VIGS]; keyed[, DOING] <- 7 - keyed[, DOING]
alpha_of <- function(cols) { Z <- keyed[, cols]; a <- suppressWarnings(psych::alpha(Z[complete.cases(Z), ], check.keys = FALSE, warnings = FALSE))
  list(alpha = a$total$raw_alpha, n = sum(complete.cases(Z)), items = length(cols),
       item_rest = as.list(setNames(a$item.stats$r.drop, cols))) }
ALPHA <- list(all_nine = alpha_of(VIGS), without_2_4_7 = alpha_of(setdiff(VIGS, AMBIG)))
# Inter-vignette correlations (keyed), for the record.
RV <- round(cor(keyed, use = "pairwise.complete.obs"), 3)

## 2. Re-scored without 2, 4 and 7 ---------
sep6 <- separation(setdiff(SENS, AMBIG), setdiff(DOING, AMBIG))
ta6 <- fifelse(is.na(sep6), NA_character_, fifelse(sep6 >= SF2B$separation_cut, "yes", "no"))
k6 <- sum(ta6 == "yes", na.rm = TRUE); n6 <- sum(!is.na(ta6)); k9 <- sum(TA == "yes", na.rm = TRUE); n9 <- sum(!is.na(TA))
xt <- table(registered = TA, without_2_4_7 = ta6)
po <- sum(diag(xt)) / sum(xt); pe <- sum(rowSums(xt) * colSums(xt)) / sum(xt)^2
RESC <- list(
  registered = list(n = n9, tells_apart = k9, share = k9 / n9, share_ci = wilson(k9, n9), separation_mean = mean(sep_all, na.rm = TRUE), separation_sd = sd(sep_all, na.rm = TRUE)),
  without_2_4_7 = list(n = n6, tells_apart = k6, share = k6 / n6, share_ci = wilson(k6, n6), separation_mean = mean(sep6, na.rm = TRUE), separation_sd = sd(sep6, na.rm = TRUE),
                       vignettes = c(setdiff(SENS, AMBIG), setdiff(DOING, AMBIG))),
  crosstab = as.data.table(xt), agreement = po, kappa = (po - pe) / (1 - pe),
  mcnemar_p = tryCatch(mcnemar.test(xt)$p.value, error = function(e) NA_real_),
  r_separations = cor(sep_all, sep6, use = "complete.obs"),
  separation_with_8 = list(registered = r_ci(sep_all, S$sf8), without_2_4_7 = r_ci(sep6, S$sf8)),
  separation_with_self = list(registered = r_ci(sep_all, V[, "vig_self"]), without_2_4_7 = r_ci(sep6, V[, "vig_self"])),
  sf8_by_group = list(registered = M[, .(n = .N, sf8 = mean(S$sf8[.I], na.rm = TRUE)), by = .(tells_apart)],
                      without_2_4_7 = data.table(ta6 = ta6, sf8 = S$sf8)[, .(n = .N, sf8 = mean(sf8, na.rm = TRUE)), by = ta6]))

## 3. Moderation: two-group models ---------
sf2b_log("17: moderation models")
target_lines <- function(tg, rlab) {
  if (tg %in% SINGLE) return(list(syn = c(paste0("L_", tg, " =~ 1*", tg), paste0(tg, " ~~ 0*", tg), paste0("L_", tg, " ~~ 1*L_", tg),
                                          paste0("S ~~ ", rlab, "*L_", tg), paste0("RI ~~ 0*L_", tg)), ord = tg))
  x <- sf2b_validity_syntax(tg, D)
  x <- x[!grepl("^S ~~", x)]
  list(syn = c(x, paste0("S ~~ ", rlab, "*", tg)), ord = score_items[[tg]])
}
build <- function(tg, mode) {
  rlab <- switch(mode, free = "c(r_no, r_yes)", equal = "c(r_eq, r_eq)", pooled = "r_pooled")
  tl <- target_lines(tg, rlab)
  syn <- c(sf2b_scale_syntax(I8), tl$syn, if (mode == "free") "dif := r_yes - r_no")
  list(syn = syn, ord = unique(c(I8, tl$ord)))
}
build_h6 <- function(mode) {
  blab <- switch(mode, free = "c(b_no, b_yes)", equal = "c(b_eq, b_eq)", pooled = "b_pooled")
  list(syn = c(sf2b_scale_syntax(I8), paste0("S ~ ", blab, "*grp"), "RI ~ 0*grp", if (mode == "free") "dif := b_yes - b_no"), ord = I8)
}
TASKS <- c(unlist(lapply(c(TARGETS, "h6"), function(tg) lapply(c("free", "equal", "pooled"), function(md) list(target = tg, mode = md))), recursive = FALSE))
run_task <- function(t) {
  m <- if (t$target == "h6") build_h6(t$mode) else build(t$target, t$mode)
  vars <- c(m$ord, if (t$target == "h6") "grp")
  d <- DAT[!is.na(DAT$.ta), c(vars, ".ta"), drop = FALSE]
  if (t$target == "h6") d <- d[!is.na(d$grp), , drop = FALSE]
  mg <- sf2b_merge_sparse(as.matrix(d[, m$ord, drop = FALSE]), d$.ta)
  d[, m$ord] <- as.data.frame(mg$X)
  d <- d[order(d$.ta), , drop = FALSE]
  grp_arg <- if (t$mode == "pooled") NULL else ".ta"
  f <- sf2b_try(sf2b_sem(m$syn, d, ordered = m$ord, group = grp_arg))
  base <- list(task = t, n = nrow(d), n_by_group = as.list(table(d$.ta)), merges = if (nrow(mg$log)) cbind(target = t$target, mode = t$mode, mg$log) else NULL)
  if (sf2b_failed(f)) return(c(base, list(ok = FALSE, error = sf2b_msg(f), improper = "")))
  gl <- if (t$mode == "pooled") "all" else lavInspect(f, "group.label")
  if (t$mode != "pooled") stopifnot(identical(gl, c("no", "yes")))
  # orientation as keyed, in every group
  L <- lavInspect(f, "std"); if (t$mode == "pooled") L <- list(L)
  scs <- if (t$target %in% c(SINGLE, "h6")) character() else t$target
  bad <- unique(unlist(lapply(L, function(Lg) {
    lam <- Lg$lambda; k8 <- ifelse(I8 %in% D$reversed, -1, 1)
    b <- if (sum(lam[I8, "S"] * k8) <= 0) "S" else character()
    for (sc in scs) { it <- score_items[[sc]]; rv <- unlist(lapply(D$validity, `[[`, "reverse"))
      if (sum(lam[it, sc] * ifelse(it %in% rv, -1, 1)) <= 0) b <- c(b, sc) }
    b })))
  pe <- parameterEstimates(f)
  est <- pe[grepl("^(r_|b_|dif)", pe$label) & !duplicated(pe$label), c("label", "est", "se", "ci.lower", "ci.upper", "pvalue")]
  c(base, list(ok = !length(bad), error = if (length(bad)) paste("not oriented as keyed:", paste(bad, collapse = ", ")) else "",
               improper = paste(improper_mg(f), collapse = "; "), est = est, fit = sf2b_fit_stats(f), fit_obj = if (t$mode %in% c("free", "equal")) f else NULL))
}
RES <- sf2b_par_lapply(TASKS, run_task, data = list(DAT = DAT, D = D, score_items = score_items, I8 = I8, SINGLE = SINGLE,
                                                     target_lines = target_lines, build = build, build_h6 = build_h6,
                                                     improper_of = improper_of, improper_list = improper_list, improper_mg = improper_mg,
                                                     not_improper = not_improper))
for (r in Filter(function(r) !r$ok, RES)) sf2b_log("FAILED: ", r$task$target, " ", r$task$mode, ": ", r$error)
for (r in Filter(function(r) r$ok && nzchar(r$improper), RES)) sf2b_log("IMPROPER: ", r$task$target, " ", r$task$mode, ": ", r$improper)
get_res <- function(tg, md) Filter(function(r) r$task$target == tg && r$task$mode == md, RES)[[1]]
FITS <- rbindlist(lapply(RES, function(r) c(list(target = r$task$target, mode = r$task$mode, n = r$n, ok = r$ok, error = r$error, improper = r$improper),
                                            if (r$ok) r$fit[c("chisq.scaled", "df", "cfi.scaled", "rmsea.scaled", "srmr")])), fill = TRUE)
MERGES <- rbindlist(lapply(RES, `[[`, "merges"), fill = TRUE)

# Observed analogues, within each group.
obs_score <- function(tg) if (tg == "vig_self") V[, "vig_self"] else S[[tg]]
fisher_diff <- function(r1, n1, r2, n2) { z <- (atanh(r1) - atanh(r2)) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3)); 2 * pnorm(-abs(z)) }
LAT04 <- fread(sf2b_out("04_latent_validity.csv"))[form == "sf8" & sample == "all" & ri_free == FALSE]
n10 <- tryCatch(fromJSON(sf2b_out("10_deviation_maia.json"))$sf8$h2_noticing$r, error = function(e) NA_real_)
MOD <- rbindlist(lapply(c(TARGETS, "h6"), function(tg) {
  fr <- get_res(tg, "free"); eq <- get_res(tg, "equal"); po <- get_res(tg, "pooled")
  pre <- if (tg == "h6") "b_" else "r_"
  g <- function(r, lab) if (isTRUE(r$ok)) { x <- r$est[r$est$label == lab, ]; if (nrow(x)) x else NULL } else NULL
  no <- g(fr, paste0(pre, "no")); ye <- g(fr, paste0(pre, "yes")); di <- g(fr, "dif"); pl <- g(po, paste0(pre, "pooled")); eqv <- g(eq, paste0(pre, "eq"))
  lrt <- if (isTRUE(fr$ok) && isTRUE(eq$ok)) tryCatch({ x <- lavTestLRT(fr$fit_obj, eq$fit_obj); list(chisq = x[2, "Chisq diff"], df = x[2, "Df diff"], p = x[2, "Pr(>Chisq)"]) },
                                                      error = function(e) list(chisq = NA_real_, df = NA_real_, p = NA_real_, error = conditionMessage(e))) else list(chisq = NA_real_, df = NA_real_, p = NA_real_)
  # observed
  if (tg == "h6") {
    ob <- lapply(c("no", "yes"), function(gv) { keep <- M$tells_apart %in% gv; y <- S$sf8; y[!keep] <- NA; sf2b_kg(y, DAT$grp) })
    dd <- data.table(y = S$sf8, g = DAT$grp, ta = as.integer(M$tells_apart == "yes"))[complete.cases(data.table(S$sf8, DAT$grp, M$tells_apart))]
    im <- summary(lm(y ~ g * ta, data = dd))$coefficients["g:ta", ]
    obs <- list(obs_no = ob[[1]]$d, obs_yes = ob[[2]]$d, obs_n_no = ob[[1]]$n_group + ob[[1]]$n_reference, obs_n_yes = ob[[2]]$n_group + ob[[2]]$n_reference,
                obs_diff_p = im[["Pr(>|t|)"]], obs_note = "observed d (07's sf2b_kg) per group; p of the practice x tells-apart interaction (OLS, raw units)")
  } else {
    y <- obs_score(tg)
    a <- r_ci(S$sf8[M$tells_apart %in% "no"], y[M$tells_apart %in% "no"]); b <- r_ci(S$sf8[M$tells_apart %in% "yes"], y[M$tells_apart %in% "yes"])
    obs <- list(obs_no = a$r, obs_yes = b$r, obs_n_no = a$n, obs_n_yes = b$n, obs_diff_p = fisher_diff(b$r, b$n, a$r, a$n),
                obs_note = "observed Pearson r with the observed 8 per group; Fisher-z test")
  }
  ref04 <- if (tg == "maia_noticing") n10 else if (tg == "h6") tryCatch(fread(sf2b_out("07_latent_difference.csv"))[form == "sf8", latent_d], error = function(e) NA_real_)
           else if (tg %in% LAT04$score) LAT04[score == tg, est] else NA_real_
  pooled <- pl$est %||% NA_real_
  se_d <- di$se %||% NA_real_
  data.table(target = tg, parameter = if (tg == "h6") "latent d, current vs never (MIMIC)" else "latent r with the 8",
             n_no = fr$n_by_group[["no"]] %||% NA, n_yes = fr$n_by_group[["yes"]] %||% NA,
             est_no = no$est %||% NA_real_, lower_no = no$ci.lower %||% NA_real_, upper_no = no$ci.upper %||% NA_real_,
             est_yes = ye$est %||% NA_real_, lower_yes = ye$ci.lower %||% NA_real_, upper_yes = ye$ci.upper %||% NA_real_,
             difference_yes_minus_no = di$est %||% NA_real_, diff_lower = di$ci.lower %||% NA_real_, diff_upper = di$ci.upper %||% NA_real_,
             diff_se = se_d, wald_p = di$pvalue %||% NA_real_,
             lrt_chisq = lrt$chisq, lrt_df = lrt$df, lrt_p = lrt$p, equal_est = eqv$est %||% NA_real_,
             detectable_difference_80 = 2.80 * se_d,
             pooled_est = pooled, reference_04_07_10 = ref04,
             observed_no = obs$obs_no, observed_yes = obs$obs_yes, observed_n_no = obs$obs_n_no, observed_n_yes = obs$obs_n_yes, observed_diff_p = obs$obs_diff_p,
             observed_note = obs$obs_note,
             detectable_observed_r_pair = if (tg != "h6" && !is.na(pooled)) sprintf("%.2f vs %.2f", tanh(atanh(pooled) - .1085), tanh(atanh(pooled) + .1085)) else NA_character_,
             free_ok = isTRUE(fr$ok), equal_ok = isTRUE(eq$ok), improper_free = fr$improper %||% "", improper_equal = eq$improper %||% "")
}), fill = TRUE)
MOD[, wald_p_holm := p.adjust(wald_p, "holm")]
z_detect <- 2.80 * sqrt(1 / (455 - 3) + 1 / (265 - 3))

## 4. Self-rating by group ---------
SELF <- rbindlist(lapply(c("no", "yes", "all"), function(gv) {
  keep <- if (gv == "all") !is.na(M$tells_apart) else M$tells_apart %in% gv
  r <- r_ci(S$sf8[keep], V[keep, "vig_self"]); rs <- sf2b_r_ci(S$sf8[keep], V[keep, "vig_self"], "spearman")
  data.table(tells_apart = gv, n = r$n, self_mean = mean(V[keep, "vig_self"], na.rm = TRUE), self_sd = sd(V[keep, "vig_self"], na.rm = TRUE),
             sf8_mean = mean(S$sf8[keep], na.rm = TRUE), sf8_sd = sd(S$sf8[keep], na.rm = TRUE),
             r = r$r, lower = r$lower, upper = r$upper, spearman = rs$r, spearman_lower = rs$lower, spearman_upper = rs$upper)
}))
SELF_TEST <- list(fisher_p_pearson = fisher_diff(SELF[tells_apart == "yes", r], SELF[tells_apart == "yes", n], SELF[tells_apart == "no", r], SELF[tells_apart == "no", n]),
                  self_mean_t_p = t.test(V[M$tells_apart %in% "yes", "vig_self"], V[M$tells_apart %in% "no", "vig_self"])$p.value,
                  latent = MOD[target == "vig_self"])

## Write ---------
fwrite(num(copy(ITEM)), sf2b_out("17_vignette_items.csv"))
fwrite(as.data.table(RV, keep.rownames = "vignette"), sf2b_out("17_vignette_correlations_keyed.csv"))
fwrite(num(copy(MOD)), sf2b_out("17_moderation.csv"))
fwrite(num(copy(FITS), 4), sf2b_out("17_moderation_fit.csv"))
if (nrow(MERGES)) fwrite(MERGES, sf2b_out("17_moderation_category_merges.csv"))
fwrite(num(copy(SELF)), sf2b_out("17_self_rating.csv"))
fwrite(RESC$crosstab, sf2b_out("17_rescore_crosstab.csv"))
out <- list(
  analysis = "17: the vignettes. EXPLORATORY, NOT REGISTERED (beyond H8 and 08's descriptives, nothing here is in osf.io/d4a5u); specified after 08's results had been seen.",
  mode = SF2B$mode, tag = SF2B$tag, written = format(Sys.time(), "%Y-%m-%d %H:%M"),
  software = list(R = R.version.string, lavaan = as.character(packageVersion("lavaan")), psych = as.character(packageVersion("psych"))),
  n = nrow(M), choices = CHOICES, vignette_kind = as.list(VK), flagged_ambiguous = AMBIG,
  items = ITEM, alpha = ALPHA, rescore = RESC, moderation = MOD, moderation_fit = FITS,
  detectable = list(observed_fisher_z_difference_80 = z_detect,
                    note = "latent: detectable_difference_80 per target in the moderation table (2.80 x SE of the difference)"),
  self_rating = list(by_group = SELF, tests = SELF_TEST),
  models = list(failed = FITS[ok == FALSE, .(target, mode, error)], improper = FITS[ok == TRUE & nzchar(improper), .(target, mode, improper)],
                category_merges = nrow(MERGES)),
  run_seconds = round(as.numeric(difftime(Sys.time(), T0, units = "secs"))))
write_json(rnd(out), sf2b_out("17_exploratory_vignettes.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
sf2b_log("17 done in ", out$run_seconds, " s; outputs ", sf2b_out("17_*"))
print(num(copy(ITEM))[, .(vignette, kind, flagged_ambiguous, mean, sd, pct_intended_side, r_rest_keyed, r_with_8, r_with_self)])
print(rnd(ALPHA[c("all_nine", "without_2_4_7")])[[1]][1:3]); print(rnd(ALPHA$without_2_4_7[1:3]))
print(rnd(RESC[c("registered", "without_2_4_7", "agreement", "kappa", "r_separations")])); print(RESC$crosstab)
print(num(copy(MOD))[, .(target, n_no, n_yes, est_no, est_yes, difference_yes_minus_no, diff_lower, diff_upper, wald_p, lrt_p, detectable_difference_80, pooled_est, reference_04_07_10)])
print(num(copy(MOD))[, .(target, observed_no, observed_yes, observed_diff_p, improper_free, improper_equal)])
print(FITS); print(num(copy(SELF))); print(rnd(SELF_TEST[1:2]))
