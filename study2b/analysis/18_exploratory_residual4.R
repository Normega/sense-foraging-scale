# 18: what is left of the 8 once all four H2 neighbours are accounted for together.
# EXPLORATORY, NOT REGISTERED.
#
# Not part of the filed registration (osf.io/d4a5u) and not in its list of exploratory
# analyses. Specified on 15 November 2026, after the confirmatory real-mode run and
# analyses 10-17 had been seen, so every number here is hypothesis-generating only. No
# registered cell, test or output is changed.
#
# Why: H3d's multiple correlation of the 8 with its four H2 neighbours (FFMQ Observing,
# SBI Savoring the Moment, BFI-2 Aesthetic Sensitivity, Brief MAIA-2 Noticing) is .898
# (04). 12 asked what is left after Observing alone; this asks what is left after all four.
#   1. One SEM: the H1 model (the 8's content factor S with its registered agreement
#      factor) with S regressed on the four neighbours as latent factors (H3d's model,
#      Noticing on its three items only, so maia_12 never enters it). S's residual
#      variance is 1. The structural part is saturated, so the model is H3d's model
#      reparameterized and must give its R.
#   2. R^2, the residual's share (1 - R^2) and the standardized regression weights, all
#      as defined parameters with delta-method CIs; the neighbours' correlations and VIFs.
#   3. Semipartial correlations of S's disturbance (the residual) with each validity
#      measure, each measure added in a model of its own, correlated with the four
#      neighbours and with the disturbance (12's construction).
#   4. Practice group (current vs never) by MIMIC, and practice frequency (0-4) as a dose.
#   5. The agreement-free variant: 1-4 with the agreement factor free to correlate with
#      the neighbours (and the added measure), still uncorrelated with S in total.
#   6. Observed-score check: OLS of the observed 8 on the four observed neighbour scores.
#   7. Reliability of the residual, from the 8's content omega and R^2.
#
# Settings as registered (Section 5): ordinal items, WLSMV, theta, missing = "pairwise",
# conditional.x = FALSE. CHOICES below are written to the JSON.
#
# Writes output/<tag>/18_*.csv and 18_exploratory_residual4.json (aggregates only).
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/18_exploratory_residual4.R

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

# The pipeline's own data: 04_validity.R's top-level expressions up to its first fit (as
# 10-17 do): A, D, M, S, DAT, score_items; nothing is fitted.
ex04 <- parse("study2b/analysis/04_validity.R", keep.source = FALSE)
first_fit <- which(vapply(ex04, function(e) grepl("^sf2b_log\\(\"fitting", deparse(e)[1]), TRUE))
stopifnot(length(first_fit) == 1)
for (i in seq_len(first_fit - 1)) eval(ex04[[i]], envir = globalenv())
stopifnot(exists("DAT"), exists("score_items"), !exists("RES"))
# 10's check for improper solutions, taken as written.
ex10 <- parse("study2b/analysis/10_deviation_maia.R", keep.source = FALSE)
is_def <- function(e, nm) is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(nm))
eval(Filter(function(e) is_def(e, "improper_of"), as.list(ex10))[[1]], envir = globalenv())
# Single-indicator factors (no residual) and lavaan's phantom factor for an observed
# covariate (x =~ x) have a standardized loading of 1 by construction: not improper (16's rule).
improper_list <- function(f) { x <- improper_of(f)$list
  x[!grepl("^std loading (L_[[:alnum:]_]+|grp) =~ ", x) & !grepl("^std loading ([[:alnum:]_]+) =~ \\1 = ", x)] }
stopifnot(identical(M$participant, S$participant), nrow(DAT) == nrow(M))

I8 <- D$sf8
NB <- SF2B_H2
stopifnot(identical(NB, c("ffmq_observing", "sbi_moment", "bfi_aesthetic", "maia_noticing")),
          identical(score_items$maia_noticing, c("maia_02", "maia_03", "maia_04")))
NBC <- c(ffmq_observing = "obs", sbi_moment = "sbi", bfi_aesthetic = "aes", maia_noticing = "not")
PAIRS <- t(combn(seq_along(NB), 2))
MAIA_BOUND <- 9   # decision 46 under theta: l^2 <= 9, std loading <= .949 (as 13 and 14)
OMEGA_8 <- fromJSON(sf2b_out("03_h1_decision.json"))$h1$omega     # the 8's content omega (agreement as error)
ALPHA_8 <- fromJSON(sf2b_out("03_h1_decision.json"))$h1$alpha

DAT$grp <- fifelse(M$practice_group == "current", 1L, fifelse(M$practice_group == "never", 0L, NA_integer_))
DAT$freq5 <- M$freq_code
DAT$separation <- M$separation
DAT$vig_self <- M$vig_self
stopifnot(sum(DAT$grp == 1, na.rm = TRUE) == 443, sum(DAT$grp == 0, na.rm = TRUE) == 164)

# Targets, one model per measure (04's measure models; the neighbours' own measures less
# the neighbour). BFI-2 Open-Mindedness total is not added: its factor would hold the four
# Aesthetic Sensitivity items already in the model as the Aesthetic factor (double counting).
TARGETS <- list(
  ffmq_other = list(scores = c("ffmq_describing", "ffmq_acting", "ffmq_nonjudging", "ffmq_nonreactivity")),
  bfi_other = list(scores = c("bfi_intellectual", "bfi_creative")),
  eq = list(scores = "eq_decentering"), samd = list(scores = "samd_md"), ncs = list(scores = "ncs6"),
  bidr = list(scores = c("bidr_sde", "bidr_im")), bsss = list(scores = "bsss8"), spsq = list(scores = "spsq_sd"),
  asi = list(scores = "asi_pc"),
  maia_other = list(scores = c("maia_notdistracting", "maia_notworrying", "maia_attnregulation", "maia_emoawareness",
                               "maia_selfregulation", "maia_bodylistening", "maia_trusting")),
  pss = list(scores = "pss4"), phq = list(scores = "phq4_total"),
  lifesat = list(single = "lifesat"), crit_pause_days = list(single = "crit_pause_days"),
  engaged_pause = list(single = "engaged_pause"), crit_notice_first = list(single = "crit_notice_first"),
  crit_outdoors = list(single = "crit_outdoors"), separation = list(cont = "separation"), vig_self = list(single = "vig_self"))
ORDER <- c(unlist(lapply(TARGETS, function(x) c(x$scores, x$single, x$cont))))
LABELS <- rbind(SF2B_PREDICTED[, .(score, label)],
                data.table(score = c("separation", "vig_self"), label = c("Vignettes: separation score", "Vignettes: self-rating")))
VARIANTS <- c("registered", "agreement_free")

CHOICES <- list(
  model = "The H1 model (the 8: content factor S with loadings free, the agreement factor RI with loadings 1 on the unrecoded statements and variance free, as registered) plus the four H2 neighbours as latent factors with 04's syntax (variance 1, loadings started by keying, Brief MAIA-2 Noticing on maia_02-04 only, so maia_12 never enters), and S ~ b_obs*Obs + b_sbi*SBI + b_aes*Aes + b_not*Noticing with S's residual variance 1. The neighbours' correlations are free and labelled. The structural part is saturated: the model is 04's H3d model reparameterized, so its R must equal H3d's.",
  r2 = "vexp = b' Phi b (Phi the neighbours' correlation matrix); R^2 = vexp / (vexp + 1); residual share = 1 / (vexp + 1); R = sqrt(R^2); standardized weight beta_k = b_k / sqrt(vexp + 1). All are lavaan defined parameters with delta-method 95% CIs. VIF_k = the k-th diagonal of Phi^-1 (point estimates).",
  residual = "The residual is S's disturbance (variance 1), uncorrelated with the four neighbours by construction. Its correlation with a measure T is the semipartial correlation (T not residualized on the neighbours).",
  semipartial_models = "Each measure in a model of its own: the model above plus the measure's factor(s) as 04 builds them (subscales of one measure as correlated factors), every factor correlated freely with the four neighbours (c_k_T, labelled) and with S's disturbance (res_T); the structural part stays saturated, so each model's fit is that of the model with all factors simply correlated. Reported: res_T with its Wald CI, and S's total correlation with T in the same model, tot_T = (sum_k b_k c_k_T + res_T) / sqrt(vexp + 1), beside 04's cell.",
  ffmq = "FFMQ Describing, Acting with awareness, Nonjudging and Nonreactivity as one block (04's FFMQ model without Observing, which is in the model as a neighbour).",
  bfi = "BFI-2 Intellectual Curiosity and Creative Imagination as one block. The Open-Mindedness total is not added: as a factor of all 12 BFI-2 items it would hold the four Aesthetic Sensitivity items, already the Aesthetic neighbour, so they would load on two factors and Aesthetic would be counted twice. Skipped, latent and observed.",
  maia = "Brief MAIA-2's seven other subscales as one block, maia_12 kept in Not-Worrying under decision 46's bound (l^2 < 9 under theta, std loading <= .949; the bound 10 used and 13-14 apply), as 14 enters it. Fallback, used only if a block fails or is improper: each subscale alone; Not-Worrying alone without maia_12 if it then fails (decision 46's step 3).",
  bidr = "BIDR SDE and IM together as correlated factors (their registered treatment).",
  single_items = "Life satisfaction, past-week pause days, notice first, time outdoors, the engaged pause and the vignette self-rating as single-indicator ordinal factors with no residual variance (04's treatment of single items; the self-rating as 17 enters it). The vignette separation score (continuous, mean sensing - mean doing) as a single-indicator continuous factor: loading free, residual variance 0, factor variance 1. Neither vignette score has a 04 cell, so for them only this model's total r is given.",
  agreement = "Registered form (primary): RI uncorrelated with S, with the neighbours, with every added factor and with practice (RI ~ 0*x). Agreement-free variant: RI free to correlate with the four neighbours and with every added factor (04's ri_free), and still uncorrelated with S in total (04 and 14's ri_free keep S ~~ 0*RI): because S is regressed on the neighbours, this is imposed by freeing the disturbance's covariance with RI and constraining it to z_ri == -(b' ri), so Cov(S, RI) = 0. RI stays uncorrelated with practice in both forms.",
  practice_group = "Current (1, n 443) vs never (0, n 164) as 07 codes them; MIMIC (observed exogenous covariate, conditional.x = FALSE), cases with a group only. S ~ b'N + g*grp, each neighbour N_k ~ a_k*grp with residual variance 1 and residual correlations free (labelled), S's disturbance uncorrelated with them. d_residual = g, in within-group SDs of the residual. In the 8's SD: g / sqrt(vexp_w + 1), vexp_w = b' Phi_w b with Phi_w the within-group neighbour correlations, i.e. the 8's within-group SD (12's convention), and also g / sqrt(vexp_w + 1 + t^2 p(1 - p)) in the 8's total SD over the 607 (t = g + b'a, the 8's whole group difference; p = 443/607). d_8_whole = t / sqrt(vexp_w + 1) is the 8's unadjusted latent d in the same model (07's MIMIC, .831 in 12).",
  dose = "Practice frequency, the registered H5b code 0-4 (01's freq_code; 719 with an answer), as a linear observed covariate in the same MIMIC in place of the group: the slope per step in the residual's within-level SD, in the 8's within-level SD, and the 8's whole slope per step (16's latent linear trend, .317, for comparison).",
  observed = "OLS of the observed 8 on the four observed neighbour scores (complete cases): R^2 and standardized weights; the residual's Spearman correlation with each measure's observed score (Fisher z CI with the Bonett-Wright SE, sf2b_r_ci), beside the observed 8's. Practice: OLS of the 8 on the group plus the four scores, the coefficient over the model's residual SD (the residual's within-group SD) and over the pooled within-group SD of the 8 (07's d); dose: the same with freq5. Observed scores carry measurement error, so partialling them under-corrects: the observed residual still holds the part of the neighbours' true scores their observed scores miss.",
  reliability = "Reliability of a residual score: if the 8's total score has reliable content variance omega (03's content omega, agreement counted as error) and the neighbours' true scores explain R^2 of the 8's content (the latent R^2), the residual's reliable variance is omega (1 - R^2) and its error variance is still 1 - omega, so rel = omega (1 - R^2) / (1 - omega R^2). CI: the formula at R^2's CI bounds (it is monotone in R^2). Sensitivity: alpha (.769) in place of omega. It is the reliability of an ideal residual that partials the neighbours' true scores; a residual on observed neighbour scores holds error from them too and is less reliable.",
  restart = "A model that does not converge, or converges to an improper solution, is refitted once with start values from the converged base model of the same variant (lavaan start = fit, matched by parameter); everything else unchanged. The refit is kept when it converges and (for a first attempt that converged) is proper; otherwise the first attempt stands. Recorded per model (restarted_from_base), with the first attempt's failure in its warnings. Found necessary in two agreement-free MIMICs: the practice-group model did not converge from lavaan's default starts, and the dose model converged to a non-positive-definite factor covariance matrix (R^2 = 1); from the base fit's starts both are proper.",
  improper = "10's improper_of() on every fitted model (std loading >= 1, negative residual or factor variances, factor correlations >= 1, non-positive-definite factor covariance); single-indicator and covariate-phantom loadings of 1 are by construction and not listed. Warnings are recorded; every MIMIC draws lavaan's vcov warning of a smallest eigenvalue near machine zero from the observed covariate under conditional.x = FALSE (12, 07), recorded and not treated as a failure.")

# A copy of sf2b_sem that does not swallow warnings (12's), so warnings are recorded.
sem_txt <- deparse(sf2b_sem)
stopifnot(sum(grepl("suppressWarnings(do.call(lavaan::sem, args))", sem_txt, fixed = TRUE)) == 1)
sf2b_sem_loud <- eval(parse(text = sub("suppressWarnings(do.call(lavaan::sem, args))", "do.call(lavaan::sem, args)", sem_txt, fixed = TRUE)))
environment(sf2b_sem_loud) <- environment(sf2b_sem)
fit_w <- function(syn, d, ord, start = NULL) {
  w <- character()
  f <- withCallingHandlers(sf2b_try(if (is.null(start)) sf2b_sem_loud(syn, d, ordered = ord) else sf2b_sem_loud(syn, d, ordered = ord, start = start)),
                           warning = function(m) { w <<- c(w, conditionMessage(m)); invokeRestart("muffleWarning") })
  list(fit = f, warnings = unique(trimws(gsub("[[:space:]]+", " ", w))))
}

## Syntax ---------
# The neighbours' correlation (or, in a MIMIC, residual correlation) between j and k.
pl <- function(j, k) { if (j == k) return("1"); ij <- sort(c(j, k)); sprintf("p_%s_%s", NBC[ij[1]], NBC[ij[2]]) }
vexp_expr <- function() paste(c(sprintf("b_%s^2", NBC), sprintf("2*b_%s*b_%s*%s", NBC[PAIRS[, 1]], NBC[PAIRS[, 2]],
                                                                  apply(PAIRS, 1, function(x) pl(x[1], x[2])))), collapse = " + ")
core_syntax <- function(variant, covariate = NULL) {
  free <- variant == "agreement_free"
  sc <- sf2b_scale_syntax(I8)
  if (free) sc <- c(sub("^S ~~ 0\\*RI$", "S ~~ z_ri*RI", sc), sprintf("z_ri == -(%s)", paste(sprintf("b_%s*ri_%s", NBC, NB), collapse = " + ")))
  stopifnot(!free || any(sc == "S ~~ z_ri*RI"))
  nb <- sf2b_validity_syntax(NB, D, ri_free = free); nb <- nb[!grepl("^S ~~", nb)]
  s <- c(sc, nb, sprintf("%s ~~ %s*%s", NB[PAIRS[, 1]], apply(PAIRS, 1, function(x) pl(x[1], x[2])), NB[PAIRS[, 2]]),
         sprintf("S ~~ 0*%s", NB),
         paste0("S ~ ", paste(sprintf("b_%s*%s", NBC, NB), collapse = " + "), if (!is.null(covariate)) paste0(" + g*", covariate)),
         paste0("vexp := ", vexp_expr()), "r2 := vexp / (vexp + 1)", "unexplained := 1 / (vexp + 1)", "R := sqrt(vexp / (vexp + 1))",
         sprintf("beta_%s := b_%s / sqrt(vexp + 1)", NBC, NBC),
         # S's correlation with each neighbour (in a MIMIC: within-group)
         vapply(seq_along(NB), function(k) sprintf("totnb_%s := (%s) / sqrt(vexp + 1)", NBC[k],
                                                   paste(vapply(seq_along(NB), function(j) if (j == k) sprintf("b_%s", NBC[j]) else sprintf("b_%s*%s", NBC[j], pl(j, k)), ""), collapse = " + ")), ""))
  if (!is.null(covariate)) s <- c(s, paste0("RI ~ 0*", covariate), sprintf("%s ~ a_%s*%s", NB, NBC, covariate))
  s
}
target_syntax <- function(tg, variant, drop12 = FALSE) {
  free <- variant == "agreement_free"; s <- character(); ord <- character(); facs <- character(); scores <- character()
  if (length(tg$scores)) {
    Dx <- D; if (drop12) Dx$validity$maia$scores$maia_notworrying <- setdiff(D$validity$maia$scores$maia_notworrying, "maia_12")
    x <- sf2b_validity_syntax(tg$scores, Dx, ri_free = free); x <- x[!grepl("^S ~~", x)]
    if ("maia_notworrying" %in% tg$scores && !drop12) x <- c(x, "maia_notworrying =~ l_maia12*maia_12", sprintf("l_maia12^2 < %s", MAIA_BOUND))
    s <- c(s, x); ord <- c(ord, unlist(sf2b_score_items(Dx)[tg$scores])); facs <- c(facs, tg$scores); scores <- tg$scores
  }
  for (v in tg$single) {
    s <- c(s, paste0("L_", v, " =~ 1*", v), paste0(v, " ~~ 0*", v), paste0("L_", v, " ~~ 1*L_", v),
           paste0("RI ~~ ", if (free) paste0("ri_", v, "*") else "0*", "L_", v))
    ord <- c(ord, v); facs <- c(facs, paste0("L_", v))
  }
  for (v in tg$cont) {
    s <- c(s, sprintf("L_%s =~ NA*%s + start(%.3f)*%s", v, v, sd(DAT[[v]], na.rm = TRUE), v), paste0(v, " ~~ 0*", v), paste0("L_", v, " ~~ 1*L_", v),
           paste0("RI ~~ ", if (free) paste0("ri_", v, "*") else "0*", "L_", v))
    facs <- c(facs, paste0("L_", v))
  }
  nm <- sub("^L_", "", facs)
  for (i in seq_along(facs))
    s <- c(s, sprintf("S ~~ res_%s*%s", nm[i], facs[i]), sprintf("%s ~~ c_%s_%s*%s", NB, NBC, nm[i], facs[i]),
           sprintf("tot_%s := (%s + res_%s) / sqrt(vexp + 1)", nm[i], paste(sprintf("b_%s*c_%s_%s", NBC, NBC, nm[i]), collapse = " + "), nm[i]))
  list(syn = s, ord = ord, facs = facs, scores = scores, D = if (drop12) Dx else D)
}
build <- function(t) {
  if (t$kind == "base") return(list(syn = core_syntax(t$variant), ord = c(I8, unlist(score_items[NB])), scores = character(), D = D, data = "all"))
  if (t$kind == "mimic") {
    syn <- c(core_syntax(t$variant, t$covariate),
             if (t$covariate == "grp") c(sprintf("tot := g + %s", paste(sprintf("b_%s*a_%s", NBC, NBC), collapse = " + ")),
                                         "d_residual := g", "d_res_in_8_sd := g / sqrt(vexp + 1)",
                                         sprintf("d_res_in_8_total_sd := g / sqrt(vexp + 1 + tot^2*%.6f)", t$pq),
                                         "d_8_whole := tot / sqrt(vexp + 1)", sprintf("d_8_whole_in_total_sd := tot / sqrt(vexp + 1 + tot^2*%.6f)", t$pq))
             else c(sprintf("tot := g + %s", paste(sprintf("b_%s*a_%s", NBC, NBC), collapse = " + ")),
                    "slope_residual := g", "slope_res_in_8_sd := g / sqrt(vexp + 1)", "slope_8_whole := tot / sqrt(vexp + 1)"))
    return(list(syn = syn, ord = c(I8, unlist(score_items[NB])), scores = character(), D = D, data = t$covariate))
  }
  ts <- target_syntax(t$tg, t$variant, isTRUE(t$drop12))
  list(syn = c(core_syntax(t$variant), ts$syn), ord = unique(c(I8, unlist(score_items[NB]), ts$ord)), scores = ts$scores, D = ts$D, data = "all", facs = ts$facs)
}
run_task <- function(t) {
  m <- build(t)
  d <- switch(m$data, all = DAT, grp = DAT[!is.na(DAT$grp), , drop = FALSE], freq5 = DAT[!is.na(DAT$freq5), , drop = FALSE])
  t0 <- Sys.time(); fw <- fit_w(m$syn, d, m$ord); f <- fw$fit; restarted <- FALSE
  # A fit that does not converge is retried once from the converged base model of the same
  # variant (start values matched by parameter), and the retry is recorded.
  first_bad <- if (sf2b_failed(f)) sf2b_msg(f) else paste(improper_list(f), collapse = "; ")
  if (nzchar(first_bad) && t$kind != "base") {
    b <- build(list(kind = "base", variant = t$variant))
    fb <- sf2b_try(sf2b_sem(b$syn, DAT, ordered = b$ord))
    if (!sf2b_failed(fb)) { fw2 <- fit_w(m$syn, d, m$ord, start = fb)
      # Kept when it converges and, if the first converged, is proper (otherwise the first stands).
      if (!sf2b_failed(fw2$fit) && (sf2b_failed(f) || !length(improper_list(fw2$fit)))) { restarted <- TRUE
        fw <- list(fit = fw2$fit, warnings = c(paste("first attempt:", first_bad), fw$warnings, fw2$warnings)); f <- fw$fit } }
  }
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  base <- list(task = t, ok = FALSE, n = nrow(d), warnings = fw$warnings, secs = secs, improper = "", restarted = restarted)
  if (sf2b_failed(f)) return(modifyList(base, list(error = sf2b_msg(f))))
  bad <- sf2b_orientation(f, I8, c(NB, m$scores), m$D)
  pt <- parTable(f)
  # S's disturbance must be uncorrelated with the neighbours (no automatic S ~~ N).
  free_sn <- pt[pt$op == "~~" & ((pt$lhs == "S" & pt$rhs %in% NB) | (pt$rhs == "S" & pt$lhs %in% NB)) & pt$free > 0, ]
  if (nrow(free_sn)) bad <- c(bad, "S's disturbance free with a neighbour")
  pe <- parameterEstimates(f); ss <- standardizedSolution(f)
  est <- pe[nzchar(pe$label) & !grepl("^l_maia12$", pe$label), c("label", "op", "est", "se", "ci.lower", "ci.upper", "pvalue")]
  ri <- ss[grepl("^(ri_|z_ri$)", ss$label), c("label", "est.std", "se", "ci.lower", "ci.upper", "pvalue")]
  extra <- list()
  std <- lavInspect(f, "std")$lambda
  for (v in t$tg$cont %||% character()) {   # the continuous single indicator's orientation
    lv <- std[v, paste0("L_", v)]; extra$cont_loading <- lv
    if (lv < 0) bad <- c(bad, paste0("L_", v))
  }
  if ("maia_notworrying" %in% m$scores && !isTRUE(t$drop12)) extra$maia12_std_loading <- std["maia_12", "maia_notworrying"]
  phi <- lavInspect(f, "cor.lv")[NB, NB]
  modifyList(base, c(list(ok = !length(bad), error = if (length(bad)) paste("not oriented as keyed:", paste(bad, collapse = ", ")) else "",
                          est = est, ri = ri, phi = phi, fit = sf2b_fit_stats(f), improper = paste(improper_list(f), collapse = "; ")), extra))
}

## Tasks ---------
G <- DAT[!is.na(DAT$grp), , drop = FALSE]
PQ <- mean(G$grp) * (1 - mean(G$grp))
TASKS <- list()
add <- function(...) TASKS[[length(TASKS) + 1L]] <<- list(...)
for (v in VARIANTS) {
  add(kind = "base", variant = v, name = "base")
  for (n in names(TARGETS)) add(kind = "target", variant = v, name = n, tg = TARGETS[[n]])
  add(kind = "mimic", variant = v, name = "practice group", covariate = "grp", pq = PQ)
  add(kind = "mimic", variant = v, name = "practice frequency (dose)", covariate = "freq5")
}
tid <- function(t) paste(t$variant, t$name, if (isTRUE(t$drop12)) "(maia_12 out)" else "", sep = " | ")
PAR <- list(DAT = DAT, D = D, score_items = score_items, I8 = I8, NB = NB, NBC = NBC, PAIRS = PAIRS, MAIA_BOUND = MAIA_BOUND,
            pl = pl, vexp_expr = vexp_expr, core_syntax = core_syntax, target_syntax = target_syntax, build = build,
            fit_w = fit_w, sf2b_sem_loud = sf2b_sem_loud, improper_of = improper_of, improper_list = improper_list)
sf2b_log("18: ", length(TASKS), " models on ", SF2B$cores, " workers")
t1 <- Sys.time()
RES <- sf2b_par_lapply(TASKS, run_task, data = PAR)
# Fallback for a block that fails or is improper: each factor alone; Not-Worrying alone
# without maia_12 if it then fails or is improper.
bad_res <- function(r) !r$ok || nzchar(r$improper)
FB <- list()
for (r in Filter(function(r) r$task$kind == "target" && bad_res(r) && length(r$task$tg$scores) > 1, RES))
  for (sc in r$task$tg$scores) FB[[length(FB) + 1]] <- list(kind = "target", variant = r$task$variant, name = paste0(r$task$name, ": ", sc),
                                                           tg = list(scores = sc), fallback_of = r$task$name)
if (length(FB)) { sf2b_log("fallback: ", length(FB), " single-factor models"); RES <- c(RES, sf2b_par_lapply(FB, run_task, data = PAR)) }
FB2 <- lapply(Filter(function(r) r$task$kind == "target" && bad_res(r) && identical(r$task$tg$scores, "maia_notworrying") && !isTRUE(r$task$drop12), RES),
              function(r) modifyList(r$task, list(drop12 = TRUE, name = paste(r$task$name, "(maia_12 out)"))))
if (length(FB2)) { sf2b_log("fallback: Not-Worrying without maia_12"); RES <- c(RES, sf2b_par_lapply(FB2, run_task, data = PAR)) }
secs_models <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
for (r in Filter(function(r) !r$ok, RES)) sf2b_log("FAILED: ", tid(r$task), ": ", r$error)
for (r in Filter(function(r) r$ok && nzchar(r$improper), RES)) sf2b_log("IMPROPER: ", tid(r$task), ": ", r$improper)

num <- function(DT, k = 3) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, function(x) round(x, k)), .SDcols = cols]; DT }
fits <- rbindlist(lapply(RES, function(r) c(list(variant = r$task$variant, kind = r$task$kind, model = r$task$name,
                                                  fallback_of = r$task$fallback_of %||% NA_character_, n = r$n, ok = r$ok, error = r$error,
                                                  improper = r$improper, restarted_from_base = r$restarted, maia12_std_loading = r$maia12_std_loading %||% NA_real_,
                                                  warnings = paste(r$warnings, collapse = " | "), seconds = round(r$secs, 1)),
                                             if (r$ok) r$fit)), fill = TRUE)
fwrite(fits, sf2b_out("18_model_fit.csv"))
get_res <- function(kind, variant, name = NULL) Filter(function(r) r$task$kind == kind && r$task$variant == variant && (is.null(name) || r$task$name == name), RES)
est_of <- function(r, labs) { e <- as.data.table(r$est); e[label %in% labs] }

## 1-2. R^2, weights, collinearity ---------
LAT04 <- fread(sf2b_out("04_latent_validity.csv"))[form == "sf8" & sample == "all" & ri_free == FALSE, .(score, r_04 = est)]
R_H3D <- fromJSON(sf2b_out("04_validity_decisions.json"))$sf8$h3d$R
REG <- rbindlist(lapply(VARIANTS, function(v) {
  r <- get_res("base", v)[[1]]; if (!r$ok) return(data.table(variant = v, error = r$error))
  e <- est_of(r, c("r2", "unexplained", "R", paste0("beta_", NBC), paste0("totnb_", NBC), paste0("b_", NBC), "vexp"))
  ph <- r$phi; vif <- diag(solve(ph))
  rbind(e[, .(variant = v, quantity = label, est, se, lower = ci.lower, upper = ci.upper, p = signif(pvalue, 3))],
        data.table(variant = v, quantity = paste0("vif_", NBC), est = unname(vif)), fill = TRUE)
}), fill = TRUE)
REG[, neighbour := names(NBC)[match(sub("^(beta|totnb|b|vif)_", "", quantity), NBC)]]
REG <- merge(REG, LAT04[, .(neighbour = score, r_04)], by = "neighbour", all.x = TRUE, sort = FALSE)
REG[!grepl("^totnb_", quantity), r_04 := NA_real_]
setcolorder(REG, c("variant", "quantity", "neighbour"))
fwrite(num(copy(REG)), sf2b_out("18_regression.csv"))
PHI <- rbindlist(lapply(VARIANTS, function(v) {
  r <- get_res("base", v)[[1]]; if (!r$ok) return(NULL)
  e <- est_of(r, apply(PAIRS, 1, function(x) pl(x[1], x[2])))
  e[, .(variant = v, pair = label, first = NB[match(sub("^p_([a-z]+)_.*", "\\1", label), NBC)], second = NB[match(sub("^p_[a-z]+_", "", label), NBC)],
        r = est, se, lower = ci.lower, upper = ci.upper)]
}))
fwrite(num(copy(PHI)), sf2b_out("18_neighbour_correlations.csv"))
R_check <- REG[variant == "registered" & quantity == "R", est]
sf2b_log(sprintf("R (registered form) %.4f; 04's H3d R %.3f", R_check, R_H3D))
R_reproduced <- isTRUE(abs(R_check - R_H3D) < .0015)
if (!R_reproduced) sf2b_log("WARNING: R does not reproduce H3d's")

## 3. Semipartial correlations ---------
TG <- rbindlist(lapply(Filter(function(r) r$task$kind == "target" && r$ok, RES), function(r) {
  e <- as.data.table(r$est)
  res <- e[grepl("^res_", label)][, .(score = sub("^res_", "", label), residual_r = est, lower = ci.lower, upper = ci.upper, p = pvalue)]
  tot <- e[grepl("^tot_", label)][, .(score = sub("^tot_", "", label), r_total = est)]
  ri <- as.data.table(r$ri)[grepl("^ri_", label) & !sub("^ri_", "", label) %in% NB][, .(score = sub("^ri_", "", label), r_agreement = est.std)]
  x <- merge(tot, res, by = "score")
  if (nrow(ri)) x <- merge(x, ri, by = "score", all.x = TRUE) else x[, r_agreement := NA_real_]
  cbind(variant = r$task$variant, model = r$task$name, fallback_of = r$task$fallback_of %||% NA_character_, x, improper = r$improper)
}), fill = TRUE)
# One row per score and variant: the block model when it is ok and proper, else its fallback.
TG[, primary := is.na(fallback_of)]
TG[, block_ok := !nzchar(improper)]
pick <- TG[, .SD[order(-(primary & block_ok), -block_ok, -primary)][1], by = .(variant, score)]
pick[, source_model := fifelse(primary, "measure block", paste0("fallback (", model, ")"))]
W <- merge(pick[variant == "registered", .(score, r_total, residual_r, lower, upper, p, model_reg = source_model, improper_reg = improper)],
           pick[variant == "agreement_free", .(score, r_total_af = r_total, residual_r_af = residual_r, lower_af = lower, upper_af = upper, p_af = p,
                                               r_agreement_af = r_agreement, model_af = source_model, improper_af = improper)], by = "score", all = TRUE)
W <- merge(W, LAT04, by = "score", all.x = TRUE)
W <- merge(W, LABELS, by = "score", all.x = TRUE)
W <- W[order(match(score, ORDER))]
setcolorder(W, c("score", "label", "r_04", "r_total", "residual_r", "lower", "upper", "p", "r_total_af", "residual_r_af", "lower_af", "upper_af", "p_af"))
W[, `:=`(p = signif(p, 3), p_af = signif(p_af, 3))]
fwrite(num(copy(W)), sf2b_out("18_residual_correlations.csv"))
fwrite(num(copy(TG)), sf2b_out("18_residual_correlations_all_models.csv"))

## 4. Practice group and dose ---------
PR <- rbindlist(lapply(VARIANTS, function(v) rbindlist(lapply(c("practice group", "practice frequency (dose)"), function(nm) {
  r <- get_res("mimic", v, nm)[[1]]
  if (!r$ok) return(data.table(variant = v, analysis = nm, error = r$error))
  e <- est_of(r, c("g", "tot", "d_residual", "d_res_in_8_sd", "d_res_in_8_total_sd", "d_8_whole", "d_8_whole_in_total_sd",
                   "slope_residual", "slope_res_in_8_sd", "slope_8_whole", "r2", "R", paste0("a_", NBC), paste0("beta_", NBC)))
  e[, .(variant = v, analysis = nm, n = r$n, quantity = label, est, se, lower = ci.lower, upper = ci.upper, p = signif(pvalue, 3),
        error = "", improper = r$improper, warnings = paste(r$warnings, collapse = " | "))]
}))), fill = TRUE)
fwrite(num(copy(PR)), sf2b_out("18_practice.csv"))

## 6. Observed-score check ---------
OBSV <- data.table(sf8 = S$sf8, S[, ..NB])
okc <- complete.cases(OBSV)
lm_o <- lm(reformulate(NB, "sf8"), data = OBSV[okc])
zs <- as.data.table(scale(OBSV[okc]))
lm_z <- lm(reformulate(NB, "sf8"), data = zs)
res_o <- rep(NA_real_, nrow(S)); res_o[okc] <- resid(lm_o)
obs_measures <- setdiff(ORDER, character())
OVAL <- c(lapply(setNames(nm = intersect(obs_measures, names(S))), function(s) S[[s]]),
          list(separation = M$separation, vig_self = M$vig_self))
OC <- rbindlist(lapply(names(OVAL), function(s) {
  a <- sf2b_r_ci(S$sf8, OVAL[[s]], "spearman"); b <- sf2b_r_ci(res_o, OVAL[[s]], "spearman")
  data.table(score = s, n = b$n, rho_8 = a$r, rho_residual = b$r, lower = b$lower, upper = b$upper, p = signif(b$p, 3))
}))
OC <- merge(OC, LABELS, by = "score", all.x = TRUE)[order(match(score, ORDER))]
setcolorder(OC, c("score", "label"))
fwrite(num(copy(OC)), sf2b_out("18_observed_correlations.csv"))
obs_reg <- data.table(term = names(coef(lm_z))[-1], std_weight = unname(coef(lm_z)[-1]),
                      lower = confint(lm_z)[-1, 1], upper = confint(lm_z)[-1, 2])
obs_r2 <- summary(lm_o)$r.squared
obs_known <- function(x, name) {
  d <- data.table(y = S$sf8, x = x, S[, ..NB]); d <- d[complete.cases(d)]
  if (name == "practice group") { n1 <- sum(d$x == 1); n0 <- sum(d$x == 0)
    sp <- sqrt(((n1 - 1) * var(d$y[d$x == 1]) + (n0 - 1) * var(d$y[d$x == 0])) / (n1 + n0 - 2)) } else sp <- sd(d$y)
  m0 <- lm(y ~ x, data = d); m1 <- lm(reformulate(c("x", NB), "y"), data = d)
  rb <- function(m, s) { cf <- coef(m)[["x"]]; ci <- confint(m)["x", ]; c(est = cf / s, lower = ci[[1]] / s, upper = ci[[2]] / s, p = summary(m)$coefficients["x", 4]) }
  rbind(data.table(analysis = name, n = nrow(d), quantity = "unadjusted, in the 8's SD", t(rb(m0, sp))),
        data.table(analysis = name, n = nrow(d), quantity = "adjusted for the four observed neighbours, in the residual's SD", t(rb(m1, sigma(m1)))),
        data.table(analysis = name, n = nrow(d), quantity = "adjusted for the four observed neighbours, in the 8's SD", t(rb(m1, sp))))
}
OK <- rbind(obs_known(DAT$grp, "practice group"), obs_known(DAT$freq5, "practice frequency (dose)"))
OK[, sd_note := fifelse(analysis == "practice group", "8's SD = pooled within-group SD (07's d); residual's SD = the adjusted model's residual SD",
                        "per step of freq5; 8's SD = its SD over the 719; residual's SD = the adjusted model's residual SD")]
fwrite(num(copy(OK)), sf2b_out("18_observed_practice.csv"))
fwrite(num(copy(rbind(data.table(term = "R2", std_weight = obs_r2), obs_reg, fill = TRUE))), sf2b_out("18_observed_regression.csv"))

## 7. Reliability of the residual ---------
rel_of <- function(om, r2) om * (1 - r2) / (1 - om * r2)
REL <- rbindlist(lapply(VARIANTS, function(v) {
  x <- REG[variant == v & quantity == "r2"]; if (!nrow(x)) return(NULL)
  rbindlist(lapply(list(c("omega (content; 03)", OMEGA_8), c("alpha (03)", ALPHA_8)), function(o) {
    om <- as.numeric(o[2])
    data.table(variant = v, reliability_of_8 = o[1], value_of_8 = om, r2 = x$est, r2_lower = x$lower, r2_upper = x$upper,
               residual_reliable_variance = om * (1 - x$est), residual_error_variance = 1 - om,
               residual_reliability = rel_of(om, x$est), lower = rel_of(om, x$upper), upper = rel_of(om, x$lower))
  }))
}))
fwrite(num(copy(REL)), sf2b_out("18_residual_reliability.csv"))

## Report ---------
gq <- function(v, q) { x <- REG[variant == v & quantity == q]; if (nrow(x)) as.list(x[, .(est, lower, upper)]) else NULL }
rnd <- function(x) if (is.data.frame(x)) x else if (is.list(x)) lapply(x, rnd) else if (is.numeric(x)) signif(x, 4) else x
out <- list(
  analysis = "18: what is left of the 8 once all four H2 neighbours are accounted for together. EXPLORATORY, NOT REGISTERED (osf.io/d4a5u does not include it); specified after the confirmatory run and analyses 10-17, with their results seen.",
  mode = SF2B$mode, tag = SF2B$tag, written = format(Sys.time(), "%Y-%m-%d %H:%M"),
  software = list(R = R.version.string, lavaan = as.character(packageVersion("lavaan")), psych = as.character(packageVersion("psych"))),
  n = nrow(DAT), n_current = sum(G$grp == 1), n_never = sum(G$grp == 0), n_freq = sum(!is.na(DAT$freq5)), choices = CHOICES,
  checks = list(R_registered_form = R_check, R_h3d_04 = R_H3D, R_reproduced = R_reproduced,
                r_8_observing_agreement_free_here = gq("agreement_free", "totnb_obs"),
                r_8_observing_agreement_free_14 = .78, omega_8 = OMEGA_8, alpha_8 = ALPHA_8),
  regression = REG, neighbour_correlations = PHI,
  residual_share = list(registered = gq("registered", "unexplained"), agreement_free = gq("agreement_free", "unexplained")),
  residual_correlations = W, practice = PR, observed = list(r2 = obs_r2, weights = obs_reg, correlations = OC, practice = OK),
  reliability = REL,
  models = list(fitted = length(RES), failed = fits[ok == FALSE, .(variant, model, error)],
                improper = fits[ok == TRUE & nzchar(improper), .(variant, model, improper)],
                fallbacks_run = length(FB) + length(FB2),
                skipped = list(bfi_openmindedness = "not added: its factor would contain the four Aesthetic Sensitivity items already in the model (double counting)")),
  run_seconds = list(models = round(secs_models), total = round(as.numeric(difftime(Sys.time(), T0, units = "secs")))), cores = SF2B$cores)
write_json(rnd(out), sf2b_out("18_exploratory_residual4.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
sf2b_log("18 done in ", out$run_seconds$total, " s; outputs ", sf2b_out("18_*"))
print(num(copy(REG))); print(num(copy(PHI)))
print(num(copy(W))[, .(score, r_04, r_total, residual_r, lower, upper, r_total_af, residual_r_af, lower_af, upper_af, model_reg, improper_reg)])
print(num(copy(PR))[, .(variant, analysis, quantity, est, lower, upper, p)])
print(num(copy(OC))); print(num(copy(OK))); print(num(copy(REL)))
