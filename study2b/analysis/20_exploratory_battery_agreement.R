# 20: the overlap and the general factor with agreement (acquiescence) modelled across the
# whole battery, not only in the 8. EXPLORATORY, NOT REGISTERED.
#
# Not part of the filed registration (osf.io/d4a5u) and not in its list of exploratory
# analyses. Specified on 16 November 2026, post hoc: after the confirmatory real-mode run
# and analyses 10-19 had been seen, so every number here is hypothesis-generating only.
# No registered cell, test or output is changed.
#
# Why: the registered models remove yes-saying from the 8 (its agreement factor) but not
# from the other measures. 04's ri_free variant shows the 8's agreement factor correlating
# with ASI-3 (.42), FFMQ Observing (.35), SAMD (.30), BSSS-8 (.31), PHQ-4 (.23), EQ (.20) and
# FFMQ Nonjudging (-.32), and 14 found that freeing it lowers r(8, Observing) from .865 to
# .78. Shared acquiescence may therefore inflate the overlap and the general factor.
#
# Method.
#   1. A common agreement factor A for each model's instruments: A loads on every item of
#      every multi-item instrument in the model, in the raw (unrecoded) direction; its
#      loadings are equal within an instrument and free between instruments (response
#      formats differ); A's variance is 1; A is uncorrelated with every content factor.
#      Single-item measures (life satisfaction, the past-week items, the engaged pause, the
#      vignette self-rating, the separation score) and practice group get no A loading.
#      A is identified by reversed items. A model in which no instrument contributes a
#      reversed item is not identified that way (it is fitted for the record only). The
#      fallback, used when a model does not converge, is improper, is not identified or has
#      no reversed items: A on the instruments with reversed items only; then the
#      registered form (A on the 8 alone) or, without the 8, no A (13's model).
#   2. Re-estimated with battery-wide A: (a) every latent correlation of the 8 (04's
#      cells, 04's measure models) and the five H3b contrasts; (b) H3d's multiple R, by the
#      same bootstrap as 04 (same seed and resamples) and by delta method in 18's
#      reparameterization; (c) 12's residual after FFMQ Observing; (d) 18's residual after
#      all four H2 neighbours, with the practice-group MIMIC; (e) 13's general factor
#      rebuilt from the agreement-corrected latent matrix by 13's method and checks.
#   3. A's size per instrument (its standardized loadings and share of item variance in
#      the instrument's model with the 8), its stability across the pairs it enters, and
#      a diagnostic with A free to correlate with the measure's content factors.
#   Sensitivity ("calibrated"): A's loading on each instrument fixed at its value from the
#   instrument's model with the 8, where the 8's reversed statements identify it, so pairs
#   without reversed items (e.g. EQ with SAMD) are corrected too; the matrix and the
#   general factor are rebuilt from it.
#
# Limitation found in the run (16 November 2026): in FFMQ-15 and Brief MAIA-2 the reversed
# items are whole subscales (FFMQ Acting with awareness and Nonjudging, plus one Describing
# item; Brief MAIA-2 Not-Distracting and Not-Worrying), so where A is identified only by
# them (pairs without the 8) it is confounded with those subscales' content: e.g. FFMQ
# Acting with SAMD .27 -> .80, Not-Worrying with the other Brief MAIA-2 subscales changing
# sign, three Heywood cases and a smallest eigenvalue of -.47 (13: -.04). The 8's own cells
# (parts a-d), where the 8's reversed statements on one factor anchor A, are not affected
# this way; part e's matrix is, and its general factor is read with that in mind.
#
# Settings as registered (Section 5): ordinal items, WLSMV, theta, missing = "pairwise",
# conditional.x = FALSE. Brief MAIA-2 with decision 46's bound on maia_12 (residual at
# least .10 of its latent-response variance; with A loading on it the bound is
# l_maia12^2 + a_maia^2 <= 9). CHOICES below are written to the JSON.
#
# Writes output/<tag>/20_*.csv and 20_exploratory_battery_agreement.json (aggregates only).
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/20_exploratory_battery_agreement.R

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
# 10-18 do): A, D, M, S, DAT, score_items, contrast_lines; nothing is fitted.
ex04 <- parse("study2b/analysis/04_validity.R", keep.source = FALSE)
first_fit <- which(vapply(ex04, function(e) grepl("^sf2b_log\\(\"fitting", deparse(e)[1]), TRUE))
stopifnot(length(first_fit) == 1)
for (i in seq_len(first_fit - 1)) eval(ex04[[i]], envir = globalenv())
stopifnot(exists("DAT"), exists("score_items"), exists("contrast_lines"), !exists("RES"))
# Definitions taken as written from 10, 12, 13 and 18 (top-level assignments, evaluated here).
is_def <- function(e, nm) is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(nm))
take <- function(file, nm) {
  ex <- Filter(function(e) is_def(e, nm), as.list(parse(file.path("study2b/analysis", file), keep.source = FALSE)))
  stopifnot(length(ex) >= 1)
  eval(ex[[1]][[3]], envir = globalenv())
}
improper_of <- take("10_deviation_maia.R", "improper_of")
# Single-indicator factors (no residual) and covariate phantoms (x =~ x) have a
# standardized loading of 1 by construction: not improper (18's rule).
improper_list <- function(f) { x <- improper_of(f)$list
  x[!grepl("^std loading (L_[[:alnum:]_]+|grp) =~ ", x) & !grepl("^std loading ([[:alnum:]_]+) =~ \\1 = ", x)] }
BLOCKS <- take("13_exploratory_general_factor.R", "BLOCKS")
EXTENSION <- take("13_exploratory_general_factor.R", "EXTENSION")
FOCUS <- take("13_exploratory_general_factor.R", "FOCUS")
BATTERY <- setdiff(unlist(BLOCKS, use.names = FALSE), EXTENSION)
TG12 <- take("12_exploratory_overlap.R", "TARGETS")
SINGLE12 <- take("12_exploratory_overlap.R", "SINGLE")
TG18 <- take("18_exploratory_residual4.R", "TARGETS")
NB <- SF2B_H2
NBC <- take("18_exploratory_residual4.R", "NBC")
PAIRS <- t(combn(seq_along(NB), 2))
pl <- take("18_exploratory_residual4.R", "pl")
vexp_expr <- take("18_exploratory_residual4.R", "vexp_expr")
environment(pl) <- environment(vexp_expr) <- globalenv()
for (nm in c("orient", "sl_fit", "sl_extend", "heywood", "ext1", "d_from"))
  assign(nm, take("13_exploratory_general_factor.R", nm), envir = globalenv())
stopifnot(identical(NB, c("ffmq_observing", "sbi_moment", "bfi_aesthetic", "maia_noticing")), identical(names(NBC), NB))

I8 <- D$sf8
N <- nrow(DAT)
DAT$grp <- fifelse(M$practice_group == "current", 1, fifelse(M$practice_group == "never", 0, NA_real_))
DAT$separation <- M$separation
DAT$vig_self <- M$vig_self
P_CUR <- mean(DAT$grp, na.rm = TRUE)
G_DAT <- DAT[!is.na(DAT$grp), , drop = FALSE]
PQ <- mean(G_DAT$grp) * (1 - mean(G_DAT$grp))
stopifnot(sum(DAT$grp == 1, na.rm = TRUE) == 443, sum(DAT$grp == 0, na.rm = TRUE) == 164)
MAIA_BOUND <- 9   # decision 46 under theta: residual >= .10 of the latent-response variance

# Instruments: the A-loading groups. Every validity score belongs to its instrument.
INST_OF <- unlist(lapply(names(D$validity), function(k) setNames(rep(k, length(D$validity[[k]]$scores)), names(D$validity[[k]]$scores))))
ALL_REV <- c(intersect(I8, D$reversed), unlist(lapply(D$validity, `[[`, "reverse")))
INST_LABEL <- c(sf8 = "The 8", ffmq15 = "FFMQ-15", maia = "Brief MAIA-2", bidr16 = "BIDR-16", phq4 = "PHQ-4", sbi = "SBI Savoring the Moment",
                bfi2 = "BFI-2 Open-Mindedness", eq = "EQ Decentering", pss4 = "PSS-4", bsss = "BSSS-8", ncs = "NCS-6",
                spsq = "SPSQ-SF Sensory Discomfort", asi = "ASI-3 Physical Concerns", samd = "SAMD Meaning-Detecting")
LABELS <- c(setNames(SF2B_PREDICTED$label, SF2B_PREDICTED$score), sf8 = "The 8", grp = "Practice group (current vs never)",
            separation = "Vignettes: separation score", vig_self = "Vignettes: self-rating")

CHOICES <- list(
  status = "Exploratory and post hoc: specified 16 November 2026 after the confirmatory run and analyses 10-19 were seen. Not registered.",
  agreement_factor = "A: one factor per model, loading on every item of every multi-item instrument in the model in the raw (unrecoded) direction; loadings equal within an instrument (one label per instrument, aq_<instrument>) and free between instruments; A's variance fixed at 1 (not the 8's loading), so models with and without the 8 are parameterized alike; in a model with the 8 alone this is the registered agreement factor reparameterized (equal loadings, variance 1 instead of loadings 1, variance free). Start values +.4. A uncorrelated with every content factor (A ~~ 0*f, and A ~ 0*grp in the MIMIC). A's sign is arbitrary: every reported A quantity is oriented so the 8's loading is positive (or, without the 8, the loading of the instrument with the most reversed items in the model).",
  instruments = "An instrument is a questionnaire: FFMQ-15, Brief MAIA-2, BIDR-16, PHQ-4, SBI, BFI-2 Open-Mindedness, EQ Decentering, PSS-4, BSSS-8, NCS-6, SPSQ-SF, ASI-3, SAMD, and the 8. Subscales of one instrument share its A loading (e.g. FFMQ Observing and FFMQ Describing items load aq_ffmq15).",
  single_items = "Life satisfaction, past-week pause days, notice first and time outdoors, the engaged pause and the vignette self-rating as single-indicator ordinal factors with no residual (04's treatment); the separation score as 18's continuous single-indicator factor; practice group as 13's phantom factor (matrix) or 18's observed covariate (MIMIC). None gets an A loading; A is uncorrelated with all of them.",
  reversed_items = "Reversed items per instrument from the definitions (scoring.reverse_items in study2b/materials, study2/materials) where they list them; for FFMQ-15, BIDR-16, PHQ-4 and EQ the definitions carry no list and config.R's manual keys are used (checked there against the definitions where present). Written to 20_reversed_items.csv. Identification of A is judged on the items in the model: FFMQ Observing alone or Brief MAIA-2 Noticing alone has no reversed items.",
  identification = "A model is 'identified by reversed items' when at least one instrument loading on A contributes a reversed item. A model without one is fitted for the record but never used: its A is identified only by departures from equal content loadings (in the dry run of this script EQ with SAMD gave r = 1.26 and PHQ-4 with ASI-3 an A loading of 1.58 on PHQ-4, i.e. A turned into a content factor).",
  usable = "A fit is used when it converges, every factor is oriented as keyed, it is proper (10's improper_of(); single-indicator and phantom loadings excepted), the parameter covariance matrix has no eigenvalue below -1e-6, every free A loading has a finite SE, and A is identified by reversed items.",
  fallback = "Per model, in order, until one is usable: full (A on every instrument); A on the instruments with reversed items in the model only (the 'A identified by reversed items' variant); then, with the 8 in the model, the registered form (A on the 8 alone: 04, 12, 13, 18's models), or without the 8 no A at all (13's pair model). A model whose full form needed a fallback is listed with the step used.",
  restart = "12-18-style regression models (parts c, d) that do not converge, or converge improperly, are refitted once with start values from the converged base model of the same form (18's rule).",
  maia = "Brief MAIA-2 with decision 46's bound on maia_12 in every model holding Not-Worrying. With A loading on maia_12 the bound keeps maia_12's standardized residual at least .10: l_maia12^2 + aq_maia^2 < 9 (theta); under calibration, l_maia12^2 < 9 - aq_maia^2.",
  a_correlations = "Part a: the 8 with each measure in 04's measure model (FFMQ five facets; BFI-2 three facets; the Open-Mindedness total as one factor of its 12 items; Brief MAIA-2 eight subscales; BIDR both scales; each other measure alone; single items as single indicators), with 04's contrast definitions where both measures are in the model and 04's two pair models for contrasts 3 and 5. Beside each: the registered value (04; Brief MAIA-2 from decision 46's bound, refitted here) and 04's ri_free value (the 8's agreement factor free with the measure's factors; Brief MAIA-2 refitted here with the bound, since 04's failed). Single items get no A loading, so their cells are the registered ones by construction.",
  h3d = "Part b: the 8 with the four H2 neighbours as standalone factors (04's H3d model) plus battery A over the 8, FFMQ (Observing's items), SBI, BFI-2 (Aesthetic's items) and Brief MAIA-2 (Noticing's items). R = sqrt(r' Phi^-1 r) from the latent correlations; CI from the same bootstrap as 04 (sf2b_boot_indices with SF2B$seed_boot, the same resamples; full A in every resample; a failed resample is left out and counted). Also the delta-method CI of 18's reparameterization.",
  residual_obs = "Part c: 12's residual model, direction '8 given Observing' only (S ~ b*Obs, S residual variance 1), with battery A, and 12's targets each in a model of its own (correlated with Obs and with S's disturbance): R^2 and the semipartial correlations.",
  residual_four = "Part d: 18's model (S regressed on the four neighbours, S's residual variance 1) with battery A; 18's targets each in a model of its own; R^2, standardized weights, semipartials, and the practice-group MIMIC (current 1 vs never 0, cases with a group; A ~ 0*grp). 18's block fallback (each factor of a failed block alone) is applied if a block is unusable in every form.",
  general_factor = "Part e: 13's latent matrix rebuilt with battery A in every cell (13's 140 models: four within-instrument and 136 pair models; within-instrument cells from the instrument's own model), then 13's method unchanged: psych::fa.parallel (minres, seed SF2B$seed) for the number of group factors (at least 3), Schmid-Leiman (minres, oblimin), the general factor oriented on the 8; checks: one factor fewer and more, one-factor minres, Jennrich-Bentler bifactor, the general factor without the 8 (the 8 extended by Dwyer's method), without BFI-2 Intellectual Curiosity, and k = 3-9. Pause days and practice group are extension variables (Dwyer). Net of G as 13: partial r, practice-group d net of G. 13's observed principal-component check is not repeated (observed scores cannot be corrected for A). 13's own matrix is passed through the same code as a reproduction check.",
  calibrated = "Sensitivity: A's loading on each instrument fixed at its estimate in that instrument's model with the 8 (part e's 8-pairs, oriented on the 8), in every model without the 8, so cells between two instruments without reversed items are corrected for A too, under the assumption that one battery-wide A with stable loadings (theta metric: item residual SD 1 in every model) holds. A cell whose calibrated model is unusable takes its primary value and is flagged.",
  a_size = "Part f: per instrument, from its model with the 8 (part e's 8-pair, used form): the raw A loading (theta) with its CI, A's standardized loadings on the instrument's items (mean, min, max) and A's mean share of item latent-response variance (mean squared standardized loading). Stability: the raw A loading over every used model the instrument enters. Diagnostic: the same 8-pair with A free to correlate with the measure's factors (S still uncorrelated with A).",
  improper = "10's improper_of() on every fitted model; warnings are recorded. Every model with practice group (phantom or MIMIC) draws lavaan's vcov warning of a smallest eigenvalue near machine zero from the observed covariate under conditional.x = FALSE (12, 13, 18): recorded, not treated as a failure (the -1e-6 threshold).")

# A copy of sf2b_sem that does not swallow warnings (12's and 18's), so warnings are recorded.
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

## The model builder ---------
# sp: has8; scores (validity factors); singles (ordinal single items); conts (continuous
# single items); grp_phantom; a_mode (full, reversed_only, registered, none, calibrated);
# a_fixed (calibrated loadings by instrument); a_free_with (factors whose correlation with
# A is free: diagnostic and 04's ri_free); s_corr (label S ~~ r_x*x, as 04); struct
# (structural lines); covariate (MIMIC); data.
groups_of <- function(sp) {
  g <- list(); if (isTRUE(sp$has8)) g$sf8 <- I8
  for (s in sp$scores %||% character()) { k <- INST_OF[[s]]; g[[k]] <- unique(c(g[[k]], score_items[[s]])) }
  g
}
a_on_of <- function(sp, g = groups_of(sp)) {
  revd <- vapply(g, function(it) any(it %in% ALL_REV), TRUE)
  on <- switch(sp$a_mode, full = names(g), reversed_only = names(g)[revd], registered = intersect("sf8", names(g)), none = character(),
               calibrated = names(g)[names(g) %in% names(sp$a_fixed) & !is.na(sp$a_fixed[names(g)])])
  list(on = on, revd = revd, identified = sp$a_mode %in% c("registered", "none", "calibrated") || !length(on) || any(revd[on]))
}
build_model <- function(sp) {
  sc <- sp$scores %||% character(); sg <- sp$singles %||% character(); ct <- sp$conts %||% character()
  scorr <- isTRUE(sp$has8) && isTRUE(sp$s_corr)
  syn <- character(); ord <- character(); facs <- character()
  if (isTRUE(sp$has8)) { syn <- c(syn, paste0("S =~ NA*", paste(I8, collapse = " + ")), "S ~~ 1*S"); ord <- I8; facs <- c(facs, sf8 = "S") }
  if (length(sc)) {
    x <- sf2b_validity_syntax(sc, D, ri = FALSE)
    if (!scorr) x <- x[!grepl("^S ~~", x)]
    syn <- c(syn, x); ord <- c(ord, unlist(score_items[sc])); facs <- c(facs, setNames(sc, sc))
  }
  for (v in sg) {
    syn <- c(syn, sprintf("L_%s =~ 1*%s", v, v), sprintf("%s ~~ 0*%s", v, v), sprintf("L_%s ~~ 1*L_%s", v, v), if (scorr) sprintf("S ~~ r_%s*L_%s", v, v))
    ord <- c(ord, v); facs <- c(facs, setNames(paste0("L_", v), v))
  }
  for (v in ct) {
    syn <- c(syn, sprintf("L_%s =~ NA*%s + start(%.3f)*%s", v, v, sd(DAT[[v]], na.rm = TRUE), v), sprintf("%s ~~ 0*%s", v, v), sprintf("L_%s ~~ 1*L_%s", v, v),
             if (scorr) sprintf("S ~~ r_%s*L_%s", v, v))
    facs <- c(facs, setNames(paste0("L_", v), v))
  }
  if (isTRUE(sp$grp_phantom)) { syn <- c(syn, "L_grp =~ 1*grp", "grp ~~ 0*grp", if (scorr) "S ~~ r_grp*L_grp"); facs <- c(facs, grp = "L_grp") }
  g <- groups_of(sp); ao <- a_on_of(sp, g); on <- ao$on
  fixed <- sp$a_mode == "calibrated"
  if (length(on)) {
    terms <- unlist(lapply(on, function(k) if (fixed) sprintf("%.6f*%s", sp$a_fixed[[k]], g[[k]]) else c(rbind(paste0("NA*", g[[k]]), paste0("start(.4)*", g[[k]]), paste0("aq_", k, "*", g[[k]])))))
    syn <- c(syn, paste0("A =~ ", paste(terms, collapse = " + ")), "A ~~ 1*A",
             vapply(unname(facs), function(f) if (f %in% (sp$a_free_with %||% character())) sprintf("A ~~ ra_%s*%s", f, f) else sprintf("A ~~ 0*%s", f), ""))
    if (!is.null(sp$covariate)) syn <- c(syn, sprintf("A ~ 0*%s", sp$covariate))
  }
  if ("maia_notworrying" %in% sc) {
    syn <- c(syn, "maia_notworrying =~ l_maia12*maia_12",
             if ("maia" %in% on && !fixed) sprintf("l_maia12^2 + aq_maia^2 < %s", MAIA_BOUND)
             else if ("maia" %in% on) sprintf("l_maia12^2 < %.6f", MAIA_BOUND - sp$a_fixed[["maia"]]^2)
             else sprintf("l_maia12^2 < %s", MAIA_BOUND))
  }
  syn <- c(syn, sp$struct)
  list(syn = syn, ord = unique(ord), facs = facs, groups = g, on = on, revd = ao$revd, identified = ao$identified, scores = sc)
}
orient_ok <- function(f, m) {
  if ("sf8" %in% names(m$facs)) return(sf2b_orientation(f, I8, m$scores, D))
  if (!length(m$scores)) return(character())
  L <- lavInspect(f, "std")$lambda; its <- sf2b_score_items(D); rv <- unlist(lapply(D$validity, `[[`, "reverse"))
  Filter(function(sc) sum(L[its[[sc]], sc] * ifelse(its[[sc]] %in% rv, -1, 1)) <= 0, m$scores)
}
run_spec <- function(sp) {
  m <- build_model(sp)
  d <- if (identical(sp$data, "grp")) G_DAT else DAT
  t0 <- Sys.time(); fw <- fit_w(m$syn, d, m$ord); f <- fw$fit; restarted <- FALSE
  first_bad <- if (sf2b_failed(f)) sf2b_msg(f) else paste(improper_list(f), collapse = "; ")
  if (nzchar(first_bad) && !is.null(sp$base)) {   # 18's restart from the base model of the same form
    mb <- build_model(sp$base); fb <- sf2b_try(sf2b_sem(mb$syn, DAT, ordered = mb$ord))
    if (!sf2b_failed(fb)) { fw2 <- fit_w(m$syn, d, m$ord, start = fb)
      if (!sf2b_failed(fw2$fit) && (sf2b_failed(f) || !length(improper_list(fw2$fit)))) { restarted <- TRUE
        fw <- list(fit = fw2$fit, warnings = c(paste("first attempt:", first_bad), fw$warnings, fw2$warnings)); f <- fw$fit } }
  }
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out <- list(spec = sp[setdiff(names(sp), "base")], ok = FALSE, usable = FALSE, error = "", improper = "", warnings = fw$warnings, secs = secs,
              restarted = restarted, a_on = m$on, identified = m$identified, n = nrow(d))
  if (sf2b_failed(f)) { out$error <- sf2b_msg(f); return(out) }
  bad <- orient_ok(f, m)
  imp <- improper_list(f)
  vc <- tryCatch(min(eigen(lavInspect(f, "vcov"), symmetric = TRUE, only.values = TRUE)$values), error = function(e) NA_real_)
  pe <- parameterEstimates(f); ss <- standardizedSolution(f)
  if (is.null(pe$label)) pe$label <- ""; if (is.null(ss$label)) ss$label <- ""   # a model with no labelled parameter has no label column
  aq <- pe[pe$op == "=~" & pe$lhs == "A", ]; aq <- aq[!duplicated(aq$label), c("label", "est", "se", "ci.lower", "ci.upper")]
  a_se_ok <- !nrow(aq) || sp$a_mode == "calibrated" || all(is.finite(aq$se))
  # A's orientation: the 8's loading positive, else that of the instrument with most reversed items.
  sgn <- 1
  if (nrow(aq) && sp$a_mode != "calibrated") {
    ref <- if ("sf8" %in% m$on) "sf8" else m$on[which.max(vapply(m$on, function(k) sum(m$groups[[k]] %in% ALL_REV), 0))]
    sgn <- if (aq$est[aq$label == paste0("aq_", ref)][1] < 0) -1 else 1
    if (sgn < 0) { lo <- aq$ci.lower; aq$est <- -aq$est; aq$ci.lower <- -aq$ci.upper; aq$ci.upper <- -lo }
  }
  if (nrow(aq)) aq$instrument <- sub("^aq_", "", aq$label)
  L <- lavInspect(f, "std")$lambda
  a_size <- if ("A" %in% colnames(L)) rbindlist(lapply(m$on, function(k) { l <- L[m$groups[[k]], "A"] * sgn
    data.table(instrument = k, items = length(l), reversed_items = sum(m$groups[[k]] %in% ALL_REV), std_mean = mean(l), std_min = min(l), std_max = max(l), share_mean = mean(l^2)) })) else NULL
  C <- lavInspect(f, "cor.lv"); fac <- m$facs[m$facs %in% rownames(C)]
  C <- C[fac, fac, drop = FALSE]; dimnames(C) <- list(names(fac), names(fac))
  est <- pe[nzchar(pe$label) & !grepl("^(aq_|ra_|l_maia12$)", pe$label) & !duplicated(pe$label), c("label", "op", "est", "se", "ci.lower", "ci.upper", "pvalue")]
  # The practice-group phantom's variance is free (it is var(grp)), so S ~~ r_grp is a covariance: the correlation is its standardized value.
  sg <- ss[ss$label == "r_grp", ]
  if (nrow(sg)) est[est$label == "r_grp", c("est", "se", "ci.lower", "ci.upper", "pvalue")] <- sg[1, c("est.std", "se", "ci.lower", "ci.upper", "pvalue")]
  ra <- ss[grepl("^ra_", ss$label), c("label", "est.std", "se", "ci.lower", "ci.upper", "pvalue")]
  if (nrow(ra)) { names(ra)[2] <- "est"; if (sgn < 0) { lo <- ra$ci.lower; ra$est <- -ra$est; ra$ci.lower <- -ra$ci.upper; ra$ci.upper <- -lo } }
  extra <- list()
  if ("maia_notworrying" %in% m$scores) extra$maia12_std_residual <- diag(lavInspect(f, "std")$theta)[["maia_12"]]
  phi <- if (all(NB %in% rownames(lavInspect(f, "cor.lv")))) lavInspect(f, "cor.lv")[NB, NB] else NULL
  vcov_ok <- !is.na(vc) && vc > -1e-6
  modifyList(out, c(list(ok = !length(bad), usable = !length(bad) && !length(imp) && vcov_ok && a_se_ok && m$identified,
                         error = if (length(bad)) paste("not oriented as keyed:", paste(bad, collapse = ", ")) else "",
                         improper = paste(imp, collapse = "; "), vcov_min_eigen = vc, vcov_ok = vcov_ok, a_se_ok = a_se_ok,
                         aq = aq, a_size = a_size, cor = C, est = est, ra = ra, phi = phi, fit = sf2b_fit_stats(f)), extra))
}
PAR <- list(DAT = DAT, G_DAT = G_DAT, D = D, score_items = score_items, I8 = I8, NB = NB, INST_OF = INST_OF, ALL_REV = ALL_REV, MAIA_BOUND = MAIA_BOUND,
            groups_of = groups_of, a_on_of = a_on_of, build_model = build_model, orient_ok = orient_ok, run_spec = run_spec,
            fit_w = fit_w, sf2b_sem_loud = sf2b_sem_loud, improper_of = improper_of, improper_list = improper_list)

# The fallback chain of a spec: full, reversed only, then registered (with the 8) or none;
# a step whose A-loading set equals an earlier step's is skipped.
chain_of <- function(sp) {
  modes <- c("full", "reversed_only", if (isTRUE(sp$has8)) "registered" else "none")
  out <- list(); seen <- list()
  for (md in modes) { s <- modifyList(sp, list(a_mode = md)); if (!is.null(s$base)) s$base$a_mode <- md; on <- sort(a_on_of(s)$on)
    if (any(vapply(seen, identical, TRUE, on))) next
    seen[[length(seen) + 1]] <- on; out[[length(out) + 1]] <- s }
  out
}
# Runs chains in waves: every chain's next step is fitted until one is usable.
run_chains <- function(chains, label) {
  pos <- rep(1L, length(chains)); done <- rep(FALSE, length(chains)); all_res <- vector("list", length(chains))
  for (i in seq_along(chains)) all_res[[i]] <- list()
  wave <- 0
  repeat {
    todo <- which(!done & pos <= lengths(chains))
    if (!length(todo)) break
    wave <- wave + 1
    sf2b_log("20 ", label, ": wave ", wave, ", ", length(todo), " models")
    R <- sf2b_par_lapply(lapply(todo, function(i) chains[[i]][[pos[i]]]), run_spec, data = PAR)
    for (j in seq_along(todo)) { i <- todo[j]; R[[j]]$step <- pos[i]; all_res[[i]][[length(all_res[[i]]) + 1]] <- R[[j]]
      if (isTRUE(R[[j]]$usable)) done[i] <- TRUE else pos[i] <- pos[i] + 1L }
  }
  lapply(all_res, function(rs) { u <- Filter(function(r) isTRUE(r$usable), rs)
    list(attempts = rs, used = if (length(u)) u[[1]] else NULL) })
}
# The form a used fit actually has: no A; A on the 8 alone (the registered form); A on the
# instruments with reversed items only; A on every instrument (full); calibrated.
form_of <- function(u) if (is.null(u)) "none usable" else if (!length(u$a_on)) "none (no agreement factor)" else
  if (identical(u$a_on, "sf8")) "registered (A on the 8 only)" else if (u$spec$a_mode == "reversed_only") "A on instruments with reversed items" else u$spec$a_mode
mode_used <- function(ch) form_of(ch$used)
num <- function(DT, k = 3) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, function(x) round(x, k)), .SDcols = cols]; DT }
get_est <- function(r, lab) { if (is.null(r) || is.null(r$est)) return(NULL); e <- as.data.table(r$est); e[label %in% lab] }

## Reversed items ---------
REVTAB <- rbindlist(c(list(data.table(instrument = "sf8", label = "The 8", items = length(I8), reversed = length(intersect(I8, D$reversed)),
                                      reversed_ids = paste(intersect(I8, D$reversed), collapse = " "), key_source = "definition (sf2b-scale.json scoring.reverse_items)")),
  lapply(names(D$validity), function(k) { v <- D$validity[[k]]; d <- sf2b_mat(paste0(v$slug, ".json"))
    data.table(instrument = k, label = INST_LABEL[[k]], items = length(v$items), reversed = length(v$reverse), reversed_ids = paste(v$reverse, collapse = " "),
               key_source = if (!is.null(d$scoring$reverse_items)) "definition (scoring.reverse_items)" else "config.R manual key (the definition lists none)") })))
SCORE_REV <- rbindlist(lapply(names(INST_OF), function(s) data.table(score = s, instrument = INST_OF[[s]], items = length(score_items[[s]]),
                                                                     reversed = sum(score_items[[s]] %in% ALL_REV))))
fwrite(REVTAB, sf2b_out("20_reversed_items.csv"))
fwrite(SCORE_REV, sf2b_out("20_reversed_items_by_score.csv"))

## Part e's matrix models (13's 140), part a's extra models, references ---------
spec_blocks <- function(bs, a_mode = "full", ...) {
  sc <- unlist(lapply(setdiff(bs, c("sf8", "lifesat", "pause", "grp")), function(b) BLOCKS[[b]]), use.names = FALSE)
  sg <- c(if ("lifesat" %in% bs) "lifesat", if ("pause" %in% bs) "crit_pause_days")
  has8 <- "sf8" %in% bs
  list(blocks = bs, has8 = has8, scores = sc, singles = sg, grp_phantom = "grp" %in% bs, a_mode = a_mode, s_corr = TRUE,
       struct = if (has8) contrast_lines(sc) else character(), ...)
}
multi <- names(BLOCKS)[lengths(BLOCKS) > 1]
JOBS <- c(lapply(multi, identity), combn(names(BLOCKS), 2, simplify = FALSE))
stopifnot(length(JOBS) == 140)
# Part a's other cells of the 8 (04's models not in 13's matrix).
EXTRA8 <- list(bfi_total = list(scores = "bfi_openmindedness"), engaged_pause = list(singles = "engaged_pause"),
               crit_notice_first = list(singles = "crit_notice_first"), crit_outdoors = list(singles = "crit_outdoors"),
               pair3 = list(scores = c("sbi_moment", "bfi_openmindedness")), pair5 = list(scores = c("ffmq_observing", "samd_md")))
spec_extra8 <- function(nm, a_mode = "full", ...) { x <- EXTRA8[[nm]]
  list(blocks = c("sf8", nm), has8 = TRUE, scores = x$scores %||% character(), singles = x$singles %||% character(), a_mode = a_mode, s_corr = TRUE,
       struct = contrast_lines(x$scores %||% character()), ...) }
fac_names <- function(sp) c(sp$scores, if (length(sp$singles)) paste0("L_", sp$singles), if (isTRUE(sp$grp_phantom)) "L_grp")

t_e <- Sys.time()
CH_M <- run_chains(lapply(JOBS, function(bs) chain_of(spec_blocks(bs))), "matrix (primary)")
CH_X <- run_chains(lapply(names(EXTRA8), function(nm) chain_of(spec_extra8(nm))), "part a extra")
# References: the registered form (A on the 8 only: 04/13's models) and 04's ri_free, for every model with the 8;
# the diagnostic (full A free with the measure's factors).
J8 <- Filter(function(bs) "sf8" %in% bs, JOBS)
REF_SPECS <- c(lapply(J8, function(bs) spec_blocks(bs, "registered", variant = "registered")),
               lapply(names(EXTRA8), function(nm) spec_extra8(nm, "registered", variant = "registered")),
               lapply(J8, function(bs) { s <- spec_blocks(bs, "registered", variant = "ri_free"); s$a_free_with <- fac_names(s); s }),
               lapply(setdiff(names(EXTRA8), c("pair3", "pair5")), function(nm) { s <- spec_extra8(nm, "registered", variant = "ri_free"); s$a_free_with <- fac_names(s); s }),
               lapply(Filter(function(bs) !any(bs %in% c("lifesat", "pause", "grp")), J8), function(bs) { s <- spec_blocks(bs, "full", variant = "a_free_diagnostic"); s$a_free_with <- fac_names(s); s }),
               list((function() { s <- spec_extra8("bfi_total", "full", variant = "a_free_diagnostic"); s$a_free_with <- fac_names(s); s })()))
sf2b_log("20 references: ", length(REF_SPECS), " models")
REF <- sf2b_par_lapply(REF_SPECS, run_spec, data = PAR)
secs_e1 <- as.numeric(difftime(Sys.time(), t_e, units = "secs"))

## A's loadings from each instrument's model with the 8 (for calibration and part f) ---------
used_of <- function(CH) lapply(CH, `[[`, "used")
U_M <- used_of(CH_M)
pair8_of <- function(b) { i <- which(vapply(JOBS, function(bs) setequal(bs, c("sf8", b)), TRUE)); CH_M[[i]] }
INST_BLOCK <- c(ffmq15 = "ffmq", maia = "maia", bidr16 = "bidr", phq4 = "phq", sbi = "sbi", bfi2 = "bfi", eq = "eq", pss4 = "pss",
                bsss = "bsss", ncs = "ncs", spsq = "spsq", asi = "asi", samd = "samd")
A_CAL <- vapply(names(INST_BLOCK), function(k) { u <- pair8_of(INST_BLOCK[[k]])$used
  if (is.null(u) || is.null(u$aq) || !nrow(u$aq) || !(k %in% u$aq$instrument)) NA_real_ else u$aq$est[u$aq$instrument == k] }, 0)
sf2b_log("calibration loadings: ", paste(names(A_CAL), round(A_CAL, 3), sep = " = ", collapse = ", "))

t_cal <- Sys.time()
JOBS_NO8 <- which(!vapply(JOBS, function(bs) "sf8" %in% bs, TRUE))
CAL_SPECS <- lapply(JOBS_NO8, function(i) { s <- spec_blocks(JOBS[[i]], "calibrated"); s$a_fixed <- A_CAL; s$variant <- "calibrated"; s })
sf2b_log("20 calibrated: ", length(CAL_SPECS), " models")
CAL <- sf2b_par_lapply(CAL_SPECS, run_spec, data = PAR)
secs_cal <- as.numeric(difftime(Sys.time(), t_cal, units = "secs"))

## Parts c and d: 12's and 18's residual models ---------
t_cd <- Sys.time()
# 12: S ~ b*Obs.
struct12_base <- c("S ~ b*ffmq_observing", "r2 := b^2 / (1 + b^2)", "unexplained := 1 / (1 + b^2)", "r_8_obs := b / sqrt(1 + b^2)")
spec12 <- function(tg_name, a_mode = "full") {
  tg <- if (is.null(tg_name)) character() else TG12[[tg_name]]
  sg <- intersect(tg, SINGLE12); sc <- setdiff(tg, SINGLE12)
  st <- struct12_base
  for (t in tg) { fn <- if (t %in% SINGLE12) paste0("L_", t) else t
    st <- c(st, sprintf("S ~~ res_%s*%s", t, fn), sprintf("ffmq_observing ~~ oth_%s*%s", t, fn), sprintf("tot_%s := (b*oth_%s + res_%s) / sqrt(1 + b^2)", t, t, t)) }
  base <- list(has8 = TRUE, scores = "ffmq_observing", a_mode = a_mode, struct = struct12_base)
  list(part = "c", name = tg_name %||% "base", has8 = TRUE, scores = c("ffmq_observing", sc), singles = sg, a_mode = a_mode, struct = st,
       base = if (!is.null(tg_name)) base else NULL)
}
# 18: S ~ the four neighbours.
core18 <- function(covariate = NULL) c(
  sprintf("%s ~~ %s*%s", NB[PAIRS[, 1]], apply(PAIRS, 1, function(x) pl(x[1], x[2])), NB[PAIRS[, 2]]),
  sprintf("S ~~ 0*%s", NB),
  paste0("S ~ ", paste(sprintf("b_%s*%s", NBC, NB), collapse = " + "), if (!is.null(covariate)) paste0(" + g*", covariate)),
  paste0("vexp := ", vexp_expr()), "r2 := vexp / (vexp + 1)", "unexplained := 1 / (vexp + 1)", "R := sqrt(vexp / (vexp + 1))",
  sprintf("beta_%s := b_%s / sqrt(vexp + 1)", NBC, NBC),
  vapply(seq_along(NB), function(k) sprintf("totnb_%s := (%s) / sqrt(vexp + 1)", NBC[k],
         paste(vapply(seq_along(NB), function(j) if (j == k) sprintf("b_%s", NBC[j]) else sprintf("b_%s*%s", NBC[j], pl(j, k)), ""), collapse = " + ")), ""),
  if (!is.null(covariate)) sprintf("%s ~ a_%s*%s", NB, NBC, covariate))
spec18 <- function(tg_name = NULL, a_mode = "full", tg = if (!is.null(tg_name)) TG18[[tg_name]] else NULL, mimic = FALSE) {
  base <- list(has8 = TRUE, scores = NB, a_mode = a_mode, struct = core18())
  if (mimic) {
    st <- c(core18("grp"), sprintf("tot := g + %s", paste(sprintf("b_%s*a_%s", NBC, NBC), collapse = " + ")),
            "d_residual := g", "d_res_in_8_sd := g / sqrt(vexp + 1)", sprintf("d_res_in_8_total_sd := g / sqrt(vexp + 1 + tot^2*%.6f)", PQ),
            "d_8_whole := tot / sqrt(vexp + 1)", sprintf("d_8_whole_in_total_sd := tot / sqrt(vexp + 1 + tot^2*%.6f)", PQ))
    return(list(part = "d", name = "practice group", has8 = TRUE, scores = NB, a_mode = a_mode, struct = st, covariate = "grp", data = "grp", base = base))
  }
  if (is.null(tg)) return(c(base, list(part = "d", name = "base")))
  facs <- c(tg$scores, if (length(tg$single)) paste0("L_", tg$single), if (length(tg$cont)) paste0("L_", tg$cont))
  st <- core18()
  for (fn in facs) { nm <- sub("^L_", "", fn)
    st <- c(st, sprintf("S ~~ res_%s*%s", nm, fn), sprintf("%s ~~ c_%s_%s*%s", NB, NBC, nm, fn),
            sprintf("tot_%s := (%s + res_%s) / sqrt(vexp + 1)", nm, paste(sprintf("b_%s*c_%s_%s", NBC, NBC, nm), collapse = " + "), nm)) }
  list(part = "d", name = tg_name, has8 = TRUE, scores = c(NB, tg$scores), singles = tg$single %||% character(), conts = tg$cont %||% character(),
       a_mode = a_mode, struct = st, base = base)
}
# H3d's correlation model (04's), for the bootstrap and as a check.
spec_h3d <- function(a_mode = "full") list(part = "b", name = "h3d", has8 = TRUE, scores = NB, a_mode = a_mode, s_corr = TRUE)

CD_SPECS <- c(list(spec12(NULL)), lapply(names(TG12), spec12), list(spec18(NULL)), lapply(names(TG18), spec18), list(spec18(mimic = TRUE)), list(spec_h3d()))
CH_CD <- run_chains(lapply(CD_SPECS, chain_of), "parts b-d")
# 18's block fallback: a multi-factor target unusable in every form, each factor alone.
fb_names <- character()
for (i in seq_along(CD_SPECS)) { sp <- CD_SPECS[[i]]
  if (identical(sp$part, "d") && !is.null(TG18[[sp$name]]) && length(TG18[[sp$name]]$scores) > 1 && is.null(CH_CD[[i]]$used)) fb_names <- c(fb_names, sp$name) }
if (length(fb_names)) {
  FB_SPECS <- unlist(lapply(fb_names, function(n) lapply(TG18[[n]]$scores, function(s) { x <- spec18(paste0(n, ": ", s), tg = list(scores = s)); x$fallback_of <- n; x })), recursive = FALSE)
  CH_FB <- run_chains(lapply(FB_SPECS, chain_of), "part d block fallback"); CD_SPECS <- c(CD_SPECS, FB_SPECS); CH_CD <- c(CH_CD, CH_FB)
}
# The registered form of the bases (checks: 12's R^2 .753, 18's R .898, 04's r).
REG_CD <- sf2b_par_lapply(list(spec12(NULL, "registered"), spec18(NULL, "registered"), spec18(mimic = TRUE, a_mode = "registered"), spec_h3d("registered")), run_spec, data = PAR)
secs_cd <- as.numeric(difftime(Sys.time(), t_cd, units = "secs"))

## H3d bootstrap (04's resamples, battery A) ---------
t_boot <- Sys.time()
B <- SF2B$boot_h3d
idx <- sf2b_boot_indices(nrow(DAT), B, SF2B$seed_boot)
m3d <- build_model(spec_h3d("full"))
sf2b_log("20 H3d bootstrap: ", B, " resamples")
boot <- unlist(sf2b_par_lapply(idx, function(ix) {
  f <- sf2b_try(sf2b_sem(m3d$syn, DAT[ix, , drop = FALSE], ordered = m3d$ord, fast = TRUE))
  if (sf2b_failed(f)) return(NA_real_)
  C <- lavInspect(f, "cor.lv"); r <- C["S", h2]; Sig <- C[h2, h2]
  v <- tryCatch(drop(t(r) %*% solve(Sig, r)), error = function(e) NA_real_)
  if (is.na(v) || v < 0) NA_real_ else sqrt(v)
}, data = list(m3d = m3d, DAT = DAT, h2 = NB)))
secs_boot <- as.numeric(difftime(Sys.time(), t_boot, units = "secs"))

## Model log ---------
fit_row <- function(r, part, model, chain_step = NA_integer_, used = NA) c(list(part = part, model = model, variant = r$spec$variant %||% r$spec$a_mode,
  a_mode = r$spec$a_mode, chain_step = chain_step, used = used, a_on = paste(r$a_on, collapse = " "), identified_by_reversed = r$identified,
  n = r$n, ok = r$ok, usable = r$usable, error = r$error, improper = r$improper, vcov_min_eigen = r$vcov_min_eigen %||% NA_real_,
  a_se_ok = r$a_se_ok %||% NA, restarted_from_base = r$restarted, maia12_std_residual = r$maia12_std_residual %||% NA_real_,
  warnings = paste(r$warnings, collapse = " | "), seconds = round(r$secs, 1)), if (!is.null(r$fit)) r$fit)
chain_rows <- function(CH, part, namer) rbindlist(lapply(seq_along(CH), function(i) rbindlist(lapply(CH[[i]]$attempts, function(r)
  fit_row(r, part, namer(r), r$step, identical(r, CH[[i]]$used))), fill = TRUE)), fill = TRUE)
mname <- function(r) paste(r$spec$blocks %||% r$spec$name, collapse = " + ")
FITS <- rbindlist(list(chain_rows(CH_M, "e (matrix)", mname), chain_rows(CH_X, "a (extra)", mname),
                       rbindlist(lapply(REF, function(r) fit_row(r, "reference", mname(r))), fill = TRUE),
                       rbindlist(lapply(CAL, function(r) fit_row(r, "e (calibrated)", mname(r))), fill = TRUE),
                       chain_rows(CH_CD, "b-d", function(r) paste(r$spec$part, r$spec$name)),
                       rbindlist(lapply(REG_CD, function(r) fit_row(r, "b-d reference", paste(r$spec$part, r$spec$name))), fill = TRUE)), fill = TRUE)
fwrite(FITS, sf2b_out("20_model_fit.csv"))
FALLBACKS <- rbind(
  rbindlist(lapply(CH_M, function(ch) data.table(part = "e (matrix)", model = mname(ch$attempts[[1]]), first_form_usable = isTRUE(ch$attempts[[1]]$usable),
                                                 first_form_identified = ch$attempts[[1]]$identified, form_used = mode_used(ch),
                                                 a_on_used = if (is.null(ch$used)) "" else paste(ch$used$a_on, collapse = " "),
                                                 why_first_failed = if (isTRUE(ch$attempts[[1]]$usable)) "" else if (!ch$attempts[[1]]$identified) "no reversed items in the model"
                                                 else paste(Filter(nzchar, c(ch$attempts[[1]]$error, ch$attempts[[1]]$improper, if (isFALSE(ch$attempts[[1]]$vcov_ok)) "vcov not positive definite",
                                                              if (isFALSE(ch$attempts[[1]]$a_se_ok)) "A loading SE not finite")), collapse = "; ")))),
  rbindlist(lapply(c(CH_X, CH_CD), function(ch) { r1 <- ch$attempts[[1]]
    data.table(part = r1$spec$part %||% "a (extra)", model = mname(r1), first_form_usable = isTRUE(r1$usable), first_form_identified = r1$identified,
               form_used = mode_used(ch), a_on_used = if (is.null(ch$used)) "" else paste(ch$used$a_on, collapse = " "),
               why_first_failed = if (isTRUE(r1$usable)) "" else if (!r1$identified) "no reversed items in the model"
               else paste(Filter(nzchar, c(r1$error, r1$improper, if (isFALSE(r1$vcov_ok)) "vcov not positive definite", if (isFALSE(r1$a_se_ok)) "A loading SE not finite")), collapse = "; ")) })))
fwrite(FALLBACKS, sf2b_out("20_fallbacks.csv"))

## (a) The 8's latent correlations ---------
LAT04 <- fread(sf2b_out("04_latent_validity.csv"))[form == "sf8" & sample == "all" & ri_free == FALSE, .(score, r_registered = est, lower_registered = ci.lower, upper_registered = ci.upper)]
V04 <- readRDS(sf2b_derived("validity"))
RIF04 <- rbindlist(lapply(Filter(function(z) z$ok && z$task$ri_free && z$task$form == "sf8" && z$task$sample == "all", V04$results),
                          function(z) as.data.table(z$r)[grepl("^r_", label), .(score = sub("^r_", "", label), r_ri_free_04 = est)]))
ref_of <- function(variant) Filter(function(r) identical(r$spec$variant, variant), REF)
ref_r <- function(variant) rbindlist(lapply(Filter(function(r) (r$usable || (r$ok && !nzchar(r$improper))) && !any(r$spec$blocks %in% c("pair3", "pair5")), ref_of(variant)), function(r) {
  e <- as.data.table(r$est)[grepl("^r_", label) & op == "~~"]; if (!nrow(e)) return(NULL)
  data.table(model = mname(r), score = sub("^r_", "", e$label), r = e$est, lower = e$ci.lower, upper = e$ci.upper) }))
REG_HERE <- ref_r("registered"); RIF_HERE <- ref_r("ri_free")
# Primary cells: each score from its measure's model (13's 8-pairs and part a's extras; not pair3/pair5).
prim8 <- rbindlist(lapply(c(Filter(function(ch) "sf8" %in% ch$attempts[[1]]$spec$blocks, CH_M), CH_X[!names(EXTRA8) %in% c("pair3", "pair5")]), function(ch) {
  u <- ch$used; if (is.null(u)) return(data.table(model = mname(ch$attempts[[1]]), form_used = "none usable"))
  e <- as.data.table(u$est)[grepl("^r_", label) & op == "~~"]
  data.table(model = mname(u), score = sub("^r_", "", e$label), r_battery_A = e$est, lower = e$ci.lower, upper = e$ci.upper, p = e$pvalue,
             form_used = form_of(u), a_on = paste(u$a_on, collapse = " ")) }), fill = TRUE)
prim8 <- prim8[!is.na(score)]
AT <- merge(prim8, LAT04, by = "score", all.x = TRUE)
AT <- merge(AT, REG_HERE[, .(score, r_registered_refit = r, lower_registered_refit = lower, upper_registered_refit = upper)], by = "score", all.x = TRUE)
AT <- merge(AT, RIF04, by = "score", all.x = TRUE)
AT <- merge(AT, RIF_HERE[, .(score, r_ri_free_refit = r)], by = "score", all.x = TRUE)
# Registered value: 04's cell; Brief MAIA-2 (04 failed) from the refit with decision 46's bound (10's model).
AT[, registered_source := fifelse(!is.na(r_registered), "04", "refit with decision 46's bound (as 10 and 13)")]
AT[is.na(r_registered), `:=`(r_registered = r_registered_refit, lower_registered = lower_registered_refit, upper_registered = upper_registered_refit)]
AT[, ri_free_source := fifelse(!is.na(r_ri_free_04), "04", "refit with decision 46's bound")]
AT[, r_ri_free := fifelse(!is.na(r_ri_free_04), r_ri_free_04, r_ri_free_refit)]
AT[, `:=`(change_from_registered = r_battery_A - r_registered, label = LABELS[score])]
ROLE <- c(setNames(rep("H2", 4), SF2B_H2), bfi_intellectual = "H3b partner", ffmq_describing = "H3b partner", bfi_openmindedness = "H3b partner",
          maia_notworrying = "H3b partner", samd_md = "H3b partner", bidr_sde = "H3c", bidr_im = "H3c", bsss8 = "H3e", ncs6 = "H3e", spsq_sd = "H3e", asi_pc = "H3e",
          lifesat = "outcome", pss4 = "outcome", phq4_total = "outcome", crit_pause_days = "criterion", crit_notice_first = "criterion",
          crit_outdoors = "criterion", engaged_pause = "criterion")
AT[, role := fifelse(score %in% names(ROLE), ROLE[score], "other")]
AT <- AT[order(match(score, SF2B_PREDICTED$score))]
refit_check <- AT[registered_source == "04", .(max_abs_diff_registered_refit_vs_04 = max(abs(r_registered_refit - r_registered), na.rm = TRUE))]
rif_check <- AT[!is.na(r_ri_free_04), .(max_abs_diff_ri_free_refit_vs_04 = max(abs(r_ri_free_refit - r_ri_free_04), na.rm = TRUE))]
A_OUT <- AT[, .(score, label, role, r_registered, lower_registered, upper_registered, r_ri_free, r_battery_A, lower, upper, p = signif(p, 3),
                change_from_registered, form_used, a_on, model, registered_source, ri_free_source)]
fwrite(num(copy(A_OUT)), sf2b_out("20_a_correlations_with_8.csv"))
# H3b contrasts.
dif_used <- function(k) { for (ch in c(CH_M, CH_X)) { u <- ch$used; if (is.null(u)) next; e <- get_est(u, paste0("dif_", k)); if (nrow(e)) return(cbind(e, model = mname(u), form_used = form_of(u))) }; NULL }
dif_ref <- function(k, variant) { for (r in ref_of(variant)) { e <- get_est(r, paste0("dif_", k)); if (!is.null(e) && nrow(e)) return(e) }; NULL }
H3B_REG <- fromJSON(sf2b_out("04_validity_decisions.json"))$sf8$h3b
C4_10 <- tryCatch(fromJSON(sf2b_out("10_deviation_maia.json"))$sf8$contrast_4, error = function(e) NULL)
H3B <- rbindlist(lapply(1:5, function(k) { h <- SF2B_H3B[k]; x <- dif_used(k); rr <- dif_ref(k, "registered")
  reg <- if (k == 4 && !is.null(C4_10)) list(d = C4_10$difference, l = C4_10$lower, u = C4_10$upper) else list(d = H3B_REG$difference[k], l = H3B_REG$lower[k], u = H3B_REG$upper[k])
  data.table(contrast = k, first = h$first, second = h$second, difference_registered = reg$d, lower_registered = reg$l, upper_registered = reg$u,
             difference_registered_refit = if (!is.null(rr)) rr$est else NA_real_,
             difference_battery_A = if (!is.null(x)) x$est else NA_real_, lower = if (!is.null(x)) x$ci.lower else NA_real_, upper = if (!is.null(x)) x$ci.upper else NA_real_,
             above_zero = if (!is.null(x)) isTRUE(x$ci.lower > 0) else NA, model = if (!is.null(x)) x$model else "", form_used = if (!is.null(x)) x$form_used else "") }))
fwrite(num(copy(H3B)), sf2b_out("20_a_h3b_contrasts.csv"))

## (b) H3d ---------
mult_r <- function(C) { r <- C["sf8", NB]; Sig <- C[NB, NB]; sqrt(drop(t(r) %*% solve(Sig, r))) }   # run_spec names the 8's factor "sf8"
cd_used <- function(part, name) { i <- which(vapply(CD_SPECS, function(s) identical(s$part, part) && identical(s$name, name), TRUE)); if (length(i)) CH_CD[[i[1]]]$used else NULL }
cd_first <- function(part, name) { i <- which(vapply(CD_SPECS, function(s) identical(s$part, part) && identical(s$name, name), TRUE)); if (length(i)) CH_CD[[i[1]]] else NULL }
h3d_u <- cd_used("b", "h3d")
R_h3d_A <- if (!is.null(h3d_u)) mult_r(h3d_u$cor) else NA_real_
R_h3d_reg_refit <- { r <- Filter(function(r) identical(r$spec$name, "h3d"), REG_CD)[[1]]; if (r$ok) mult_r(r$cor) else NA_real_ }
H3D_REG <- fromJSON(sf2b_out("04_validity_decisions.json"))$sf8$h3d
REG18 <- fread(sf2b_out("18_regression.csv"))
b18 <- cd_used("d", "base")
d18 <- if (!is.null(b18)) as.data.table(b18$est) else NULL
boot_ci <- unname(quantile(boot, c(.025, .975), na.rm = TRUE))
H3D <- data.table(quantity = c("R registered (04, bootstrap CI)", "R agreement free (18's agreement_free form, delta CI)", "R battery A (bootstrap CI, 04's resamples)",
                               "R battery A (18's reparameterization, delta CI)", "R registered, refitted here (check)"),
                  est = c(H3D_REG$R, REG18[variant == "agreement_free" & quantity == "R", est], R_h3d_A, if (!is.null(d18)) d18[label == "R", est] else NA, R_h3d_reg_refit),
                  lower = c(H3D_REG$lower, REG18[variant == "agreement_free" & quantity == "R", lower], boot_ci[1], if (!is.null(d18)) d18[label == "R", ci.lower] else NA, NA),
                  upper = c(H3D_REG$upper, REG18[variant == "agreement_free" & quantity == "R", upper], boot_ci[2], if (!is.null(d18)) d18[label == "R", ci.upper] else NA, NA),
                  resamples_ok = c(H3D_REG$resamples_ok, NA, sum(!is.na(boot)), NA, NA), resamples = c(H3D_REG$resamples, NA, B, NA, NA),
                  upper_below_80 = c(isTRUE(H3D_REG$upper < .80), NA, isTRUE(boot_ci[2] < .80), NA, NA),
                  form_used = c("registered", "agreement_free", form_of(h3d_u), form_of(b18), "registered"))
fwrite(num(copy(H3D)), sf2b_out("20_b_h3d.csv"))
fwrite(data.table(resample = seq_along(boot), R = boot), sf2b_out("20_b_h3d_bootstrap.csv"))

## (c) 12's residual after Observing ---------
R12_SHARE <- fread(sf2b_out("12_a_variance_share.csv"))[direction == "8|obs"]
R12 <- fread(sf2b_out("12_a_residual_correlations.csv"))[focal == "the 8, residual after FFMQ Observing"]
c_base <- cd_used("c", "base")
C_SHARE <- rbindlist(lapply(c("b", "r2", "unexplained", "r_8_obs"), function(q) { e <- get_est(c_base, q); x <- R12_SHARE[label == q]
  data.table(quantity = q, est_registered = x$est, lower_registered = x$lower, upper_registered = x$upper,
             est_battery_A = if (!is.null(e) && nrow(e)) e$est else NA_real_, lower = if (!is.null(e) && nrow(e)) e$ci.lower else NA_real_,
             upper = if (!is.null(e) && nrow(e)) e$ci.upper else NA_real_, form_used = form_of(c_base)) }))
C_RES <- rbindlist(lapply(names(TG12), function(n) { u <- cd_used("c", n)
  rbindlist(lapply(TG12[[n]], function(t) { x <- R12[score == t]
    e <- get_est(u, paste0("res_", t)); tt <- get_est(u, paste0("tot_", t))
    data.table(score = t, label = LABELS[[t]], r_total_registered = x$r_total, residual_r_registered = x$residual_r, lower_registered = x$lower, upper_registered = x$upper,
               r_total_battery_A = if (!is.null(tt) && nrow(tt)) tt$est else NA_real_, residual_r_battery_A = if (!is.null(e) && nrow(e)) e$est else NA_real_,
               lower = if (!is.null(e) && nrow(e)) e$ci.lower else NA_real_, upper = if (!is.null(e) && nrow(e)) e$ci.upper else NA_real_,
               p = if (!is.null(e) && nrow(e)) signif(e$pvalue, 3) else NA_real_, form_used = form_of(u),
               a_on = if (!is.null(u)) paste(u$a_on, collapse = " ") else "") })) }))
fwrite(num(copy(C_SHARE)), sf2b_out("20_c_variance_share.csv"))
fwrite(num(copy(C_RES)), sf2b_out("20_c_residual_correlations.csv"))

## (d) 18's residual after the four neighbours ---------
D_REG <- rbindlist(lapply(c("r2", "unexplained", "R", paste0("beta_", NBC), paste0("totnb_", NBC)), function(q) {
  e <- get_est(b18, q); x <- REG18[variant == "registered" & quantity == q]; xf <- REG18[variant == "agreement_free" & quantity == q]
  data.table(quantity = q, est_registered = x$est, lower_registered = x$lower, upper_registered = x$upper, est_agreement_free_18 = if (nrow(xf)) xf$est else NA_real_,
             est_battery_A = if (!is.null(e) && nrow(e)) e$est else NA_real_, lower = if (!is.null(e) && nrow(e)) e$ci.lower else NA_real_,
             upper = if (!is.null(e) && nrow(e)) e$ci.upper else NA_real_) }))
fwrite(num(copy(D_REG)), sf2b_out("20_d_regression.csv"))
R18 <- fread(sf2b_out("18_residual_correlations.csv"))
d_target_used <- function(t) {   # the block model if usable, else its single-factor fallback
  for (i in seq_along(CD_SPECS)) { s <- CD_SPECS[[i]]; if (!identical(s$part, "d") || identical(s$name, "base") || identical(s$name, "practice group")) next
    u <- CH_CD[[i]]$used; if (is.null(u)) next
    if (!is.null(get_est(u, paste0("res_", t))) && nrow(get_est(u, paste0("res_", t)))) return(u) }
  NULL }
ORDER18 <- unlist(lapply(TG18, function(x) c(x$scores, x$single, x$cont)), use.names = FALSE)
D_RES <- rbindlist(lapply(ORDER18, function(t) { u <- d_target_used(t); x <- R18[score == t]
  e <- get_est(u, paste0("res_", t)); tt <- get_est(u, paste0("tot_", t))
  data.table(score = t, label = LABELS[[t]] %||% t, r_total_registered = x$r_total, residual_r_registered = x$residual_r, lower_registered = x$lower, upper_registered = x$upper,
             residual_r_agreement_free_18 = x$residual_r_af,
             r_total_battery_A = if (!is.null(tt) && nrow(tt)) tt$est else NA_real_, residual_r_battery_A = if (!is.null(e) && nrow(e)) e$est else NA_real_,
             lower = if (!is.null(e) && nrow(e)) e$ci.lower else NA_real_, upper = if (!is.null(e) && nrow(e)) e$ci.upper else NA_real_,
             p = if (!is.null(e) && nrow(e)) signif(e$pvalue, 3) else NA_real_, form_used = form_of(u),
             model = if (!is.null(u)) u$spec$name else "", a_on = if (!is.null(u)) paste(u$a_on, collapse = " ") else "") }))
fwrite(num(copy(D_RES)), sf2b_out("20_d_residual_correlations.csv"))
PR18 <- fread(sf2b_out("18_practice.csv"))[analysis == "practice group"]
mim <- cd_used("d", "practice group")
D_PR <- rbindlist(lapply(c("d_residual", "d_res_in_8_sd", "d_res_in_8_total_sd", "d_8_whole", "d_8_whole_in_total_sd", "r2", "R", "g", "tot"), function(q) {
  e <- get_est(mim, q); x <- PR18[variant == "registered" & quantity == q]; xf <- PR18[variant == "agreement_free" & quantity == q]
  data.table(quantity = q, n = if (!is.null(mim)) mim$n else NA_integer_, est_registered = x$est, lower_registered = x$lower, upper_registered = x$upper,
             est_agreement_free_18 = if (nrow(xf)) xf$est else NA_real_,
             est_battery_A = if (!is.null(e) && nrow(e)) e$est else NA_real_, lower = if (!is.null(e) && nrow(e)) e$ci.lower else NA_real_,
             upper = if (!is.null(e) && nrow(e)) e$ci.upper else NA_real_, p = if (!is.null(e) && nrow(e)) signif(e$pvalue, 3) else NA_real_,
             form_used = form_of(mim)) }))
fwrite(num(copy(D_PR)), sf2b_out("20_d_practice.csv"))
chk_reg <- list(r2_12_registered_refit = get_est(REG_CD[[1]], "r2")$est, r2_12 = R12_SHARE[label == "r2", est],
                R_18_registered_refit = get_est(REG_CD[[2]], "R")$est, R_18 = REG18[variant == "registered" & quantity == "R", est],
                d_residual_18_registered_refit = get_est(REG_CD[[3]], "d_residual")$est, d_residual_18 = PR18[variant == "registered" & quantity == "d_residual", est],
                R_h3d_registered_refit = R_h3d_reg_refit, R_h3d_04 = H3D_REG$R)

## (e) The general factor ---------
VARS <- unlist(BLOCKS, use.names = FALSE)
blk <- setNames(rep(names(BLOCKS), lengths(BLOCKS)), VARS)
assemble <- function(results) {   # 13's assembly: within-instrument cells from the instrument's own model
  RL <- matrix(NA_real_, length(VARS), length(VARS), dimnames = list(VARS, VARS)); diag(RL) <- 1
  for (r in Filter(Negate(is.null), results)) {
    C <- r$cor; v <- rownames(C)
    for (i in v) for (j in v) if (i != j) { if (blk[[i]] == blk[[j]] && length(r$spec$blocks) == 2) next; RL[i, j] <- C[i, j] }
  }
  RL
}
RL_A <- assemble(U_M)
# Calibrated: 8-pairs as primary; others from the calibrated fit when usable, else the primary.
CAL_USED <- vector("list", length(JOBS)); CAL_SRC <- rep("primary (8-pair)", length(JOBS))
for (i in seq_along(JOBS)) CAL_USED[i] <- list(U_M[[i]])
for (j in seq_along(JOBS_NO8)) { i <- JOBS_NO8[j]
  if (isTRUE(CAL[[j]]$usable)) { CAL_USED[[i]] <- CAL[[j]]; CAL_SRC[i] <- "calibrated" } else CAL_SRC[i] <- "primary (calibrated unusable)" }
RL_CAL <- assemble(CAL_USED)
LM13 <- fread(sf2b_out("13_latent_matrix.csv")); RL13 <- as.matrix(LM13[, -1, with = FALSE]); rownames(RL13) <- LM13$var
stopifnot(identical(rownames(RL13), VARS))
fwrite(data.table(var = VARS, round(RL_A, 4)), sf2b_out("20_e_latent_matrix.csv"))
fwrite(data.table(var = VARS, round(RL_CAL, 4)), sf2b_out("20_e_latent_matrix_calibrated.csv"))
CELLS <- rbindlist(lapply(seq_along(JOBS), function(i) data.table(model = paste(JOBS[[i]], collapse = " + "), form_used = mode_used(CH_M[[i]]),
  a_on = if (is.null(U_M[[i]])) "" else paste(U_M[[i]]$a_on, collapse = " "), calibrated_source = CAL_SRC[i])))
fwrite(CELLS, sf2b_out("20_e_cell_sources.csv"))
# Cell-by-cell change from 13.
ut <- which(upper.tri(RL_A), arr.ind = TRUE)
CELLDIFF <- data.table(a = VARS[ut[, 1]], b = VARS[ut[, 2]], r_13 = RL13[ut], r_battery_A = RL_A[ut], r_calibrated = RL_CAL[ut])
CELLDIFF[, `:=`(change = r_battery_A - r_13, change_calibrated = r_calibrated - r_13)]
fwrite(num(copy(CELLDIFF)), sf2b_out("20_e_cell_changes.csv"))

general_factor <- function(RL, tag) {   # 13's steps 2-3, unchanged
  dropped <- VARS[apply(RL, 1, function(x) any(is.na(x)))]
  BAT <- setdiff(BATTERY, dropped); EXT <- setdiff(EXTENSION, dropped)
  Rm <- RL[BAT, BAT]; Rm <- (Rm + t(Rm)) / 2
  ev <- eigen(Rm, symmetric = TRUE, only.values = TRUE)$values
  smoothed <- min(ev) <= 0
  if (smoothed) Rm <- suppressWarnings(cor.smooth(Rm))
  set.seed(SF2B$seed)
  pa <- suppressWarnings(fa.parallel(Rm, n.obs = N, fm = "minres", fa = "fa", plot = FALSE))
  K <- max(3L, as.integer(pa$nfact))
  SL <- suppressWarnings(sl_fit(Rm, K)); SLm <- suppressWarnings(sl_fit(Rm, K - 1L)); SLp <- suppressWarnings(sl_fit(Rm, K + 1L))
  one <- suppressWarnings(fa(Rm, nfactors = 1, fm = "minres", n.obs = N))
  g_one <- orient(setNames(as.numeric(one$loadings[, 1]), rownames(Rm)))
  jb <- tryCatch(suppressWarnings(fa(Rm, nfactors = K + 1, fm = "minres", rotate = "bifactor", n.obs = N)), error = function(e) e)
  g_jb <- if (inherits(jb, "error")) setNames(rep(NA_real_, nrow(Rm)), rownames(Rm)) else {
    Lj <- unclass(jb$loadings); gc <- which.max(colSums(Lj^2)); orient(setNames(Lj[, gc], rownames(Rm))) }
  B28 <- setdiff(BAT, "sf8")
  SL28 <- suppressWarnings(sl_fit(Rm[B28, B28], K))
  if (cor(SL28$g, SL$g[B28]) < 0) { SL28$g <- -SL28$g; SL28$flip <- -SL28$flip }
  g8_ext <- unname(sl_extend(SL28, Rm["sf8", B28, drop = FALSE]))[1]
  B_noint <- setdiff(BAT, "bfi_intellectual")
  SLni <- suppressWarnings(sl_fit(Rm[B_noint, B_noint], K))
  KSWEEP <- rbindlist(lapply(3:9, function(k) { x <- suppressWarnings(sl_fit(Rm, k))
    data.table(group_factors = k, max_communality = max(x$s$sl[, "h2"]), heywood = paste(heywood(x), collapse = " "),
               rms = x$s$rms, g_share_of_common = sum(x$g^2) / sum(x$s$sl[, "h2"]),
               g_sf8 = x$g[["sf8"]], g_ffmq_observing = x$g[["ffmq_observing"]], g_eq_decentering = x$g[["eq_decentering"]],
               g_ncs6 = x$g[["ncs6"]], g_bidr_sde = x$g[["bidr_sde"]], g_bidr_im = x$g[["bidr_im"]]) }))
  LOAD <- data.table(measure = BAT, label = LABELS[BAT], sl_g = SL$g[BAT], sl_g_one_fewer = SLm$g[BAT], sl_g_one_more = SLp$g[BAT],
                     one_factor = g_one[BAT], jennrich_bentler_g = g_jb[BAT],
                     sl_g_without_8 = c(SL28$g, sf8 = g8_ext)[BAT], sl_g_without_bfi_intellectual = SLni$g[BAT],
                     sl_h2 = SL$s$sl[BAT, "h2"], sl_g_share_of_common = SL$s$sl[BAT, "p2"])
  GDEF <- list(sl = c(SL$g, if (length(EXT)) sl_extend(SL, RL[EXT, BAT, drop = FALSE])),
               sl_without_8 = c(SL28$g, sf8 = g8_ext, if (length(EXT)) sl_extend(SL28, RL[EXT, B28, drop = FALSE])),
               one_factor = c(g_one, if (length(EXT)) ext1(RL[EXT, BAT, drop = FALSE], g_one)))
  net <- rbindlist(lapply(names(GDEF), function(gd) { la <- GDEF[[gd]]
    rbindlist(lapply(setdiff(c(BAT, EXT), "sf8"), function(x) { r <- RL["sf8", x]; l8 <- la[["sf8"]]; lx <- la[[x]]
      data.table(g_definition = gd, measure = x, label = LABELS[[x]], r_with_8 = r, g_loading_8 = l8, g_loading = lx, accounted_by_G = l8 * lx,
                 share_accounted = (l8 * lx) / r, partial_r_net_of_G = (r - l8 * lx) / sqrt((1 - l8^2) * (1 - lx^2)), extension = x %in% EXT) })) }))
  grp_net <- if ("grp" %in% EXT) list(d_raw = d_from(RL["sf8", "grp"])[["d_residual_sd"]], r_8_grp = RL["sf8", "grp"], p_current = P_CUR,
    net = rbindlist(lapply(names(GDEF), function(gd) { dd <- d_from(RL["sf8", "grp"], GDEF[[gd]][["sf8"]], GDEF[[gd]][["grp"]])
      data.table(g_definition = gd, g_loading_8 = GDEF[[gd]][["sf8"]], g_loading_grp = GDEF[[gd]][["grp"]], d_net_of_G = dd[["d_residual_sd"]], d_net_of_G_in_8_within_sd = dd[["d_8_within_sd"]]) }))) else NULL
  om <- tryCatch(suppressWarnings(omega(Rm, nfactors = K, fm = "minres", n.obs = N, flip = FALSE, plot = FALSE)), error = function(e) NULL)
  list(tag = tag, dropped = dropped, battery = BAT, smallest_eigenvalue = min(ev), smoothed = smoothed, pa_nfact = pa$nfact, K = K,
       fa_eigenvalues = round(pa$fa.values, 3), LOAD = LOAD, net = net, grp_net = grp_net, KSWEEP = KSWEEP,
       g_share_of_common = sum(SL$g^2) / sum(SL$s$sl[, "h2"]), g_share_of_total = mean(SL$g^2), second_order = as.numeric(unclass(SL$s$gloading)),
       rms = SL$s$rms, omega_h = if (!is.null(om)) om$omega_h else NA_real_, heywood = heywood(SL), heywood_without_bfi_intellectual = heywood(SLni),
       one_factor_share = mean(g_one^2), jb_error = if (inherits(jb, "error")) conditionMessage(jb) else "",
       agree = round(cor(LOAD[, .(sl_g, sl_g_one_fewer, sl_g_one_more, one_factor, jennrich_bentler_g, sl_g_without_8, sl_g_without_bfi_intellectual)], use = "pairwise"), 3))
}
GF13 <- general_factor(RL13, "13 (registered agreement, reproduction check)")
L13 <- fread(sf2b_out("13_general_loadings.csv"))
repro13 <- max(abs(GF13$LOAD$sl_g[match(L13$measure, GF13$LOAD$measure)] - L13$sl_g))
sf2b_log("13 reproduced through 20's code: max |sl_g difference| ", signif(repro13, 3))
GFA <- general_factor(RL_A, "battery A")
GFC <- general_factor(RL_CAL, "battery A, calibrated")
E_LOAD <- Reduce(function(a, b) merge(a, b, by = "measure", all = TRUE), list(
  GF13$LOAD[, .(measure, label, sl_g_13 = sl_g, h2_13 = sl_h2)],
  GFA$LOAD[, .(measure, sl_g_battery_A = sl_g, sl_g_one_fewer, sl_g_one_more, one_factor, jennrich_bentler_g, sl_g_without_8, sl_g_without_bfi_intellectual, h2_battery_A = sl_h2)],
  GFC$LOAD[, .(measure, sl_g_calibrated = sl_g, h2_calibrated = sl_h2)]))
E_LOAD[, `:=`(change = sl_g_battery_A - sl_g_13, change_calibrated = sl_g_calibrated - sl_g_13)]
E_LOAD <- E_LOAD[order(-sl_g_battery_A)]
fwrite(num(copy(E_LOAD)), sf2b_out("20_e_general_loadings.csv"))
E_NET <- rbind(cbind(matrix_source = "13", GF13$net), cbind(matrix_source = "battery A", GFA$net), cbind(matrix_source = "battery A, calibrated", GFC$net))
fwrite(num(copy(E_NET)), sf2b_out("20_e_net_of_general.csv"))
E_GRP <- rbind(cbind(matrix_source = "13", d_raw = GF13$grp_net$d_raw, GF13$grp_net$net), cbind(matrix_source = "battery A", d_raw = GFA$grp_net$d_raw, GFA$grp_net$net),
               cbind(matrix_source = "battery A, calibrated", d_raw = GFC$grp_net$d_raw, GFC$grp_net$net))
fwrite(num(copy(E_GRP)), sf2b_out("20_e_practice_group_net.csv"))
fwrite(num(rbind(cbind(matrix_source = "13", GF13$KSWEEP), cbind(matrix_source = "battery A", GFA$KSWEEP), cbind(matrix_source = "battery A, calibrated", GFC$KSWEEP))),
       sf2b_out("20_e_k_sensitivity.csv"))
E_SUM <- rbindlist(lapply(list(GF13, GFA, GFC), function(g) data.table(matrix_source = g$tag, parallel_analysis_nfact = g$pa_nfact, group_factors = g$K,
  smallest_eigenvalue = g$smallest_eigenvalue, smoothed = g$smoothed, g_share_of_common = g$g_share_of_common, g_share_of_total = g$g_share_of_total,
  omega_h = g$omega_h, one_factor_share = g$one_factor_share, sl_rms = g$rms, g_sf8 = g$LOAD[measure == "sf8", sl_g],
  heywood = paste(g$heywood, collapse = " "), dropped = paste(g$dropped, collapse = " "))))
fwrite(num(copy(E_SUM)), sf2b_out("20_e_summary.csv"))

## (f) A's size per instrument ---------
F_SIZE <- rbindlist(lapply(names(INST_BLOCK), function(k) { ch <- pair8_of(INST_BLOCK[[k]]); u <- ch$used
  base <- data.table(instrument = k, label = INST_LABEL[[k]], model = paste(c("sf8", INST_BLOCK[[k]]), collapse = " + "), form_used = mode_used(ch),
                     items = REVTAB[instrument == k, items], reversed_items = REVTAB[instrument == k, reversed],
                     response_range = paste(D$validity[[k]]$range, collapse = "-"))
  if (is.null(u) || is.null(u$a_size) || !(k %in% u$a_size$instrument)) return(base)
  s <- u$a_size[instrument == k]; s8 <- u$a_size[instrument == "sf8"]; a <- u$aq[u$aq$instrument == k, ]
  cbind(base, a_raw = a$est, a_lower = a$ci.lower, a_upper = a$ci.upper, a_std_mean = s$std_mean, a_std_min = s$std_min, a_std_max = s$std_max,
        a_share_mean = s$share_mean, a8_raw = u$aq$est[u$aq$instrument == "sf8"], a8_std_mean = s8$std_mean, a8_share_mean = s8$share_mean) }), fill = TRUE)
# Stability: the instrument's raw A loading over every used model (matrix, part a, parts b-d) it enters.
all_used <- c(U_M, lapply(CH_X, `[[`, "used"), lapply(CH_CD, `[[`, "used"))
STAB <- rbindlist(lapply(Filter(function(u) !is.null(u) && !is.null(u$aq) && nrow(u$aq), all_used), function(u)
  data.table(model = paste(u$spec$blocks %||% c(u$spec$part, u$spec$name), collapse = " + "), has8 = isTRUE(u$spec$has8), as.data.table(u$aq)[, .(instrument, est)])))
STAB_SUM <- STAB[, .(models = .N, median = median(est), min = min(est), max = max(est), models_without_8 = sum(!has8),
                     median_without_8 = if (any(!has8)) median(est[!has8]) else NA_real_), by = instrument]
F_SIZE <- merge(F_SIZE, STAB_SUM[, .(instrument, stability_models = models, a_raw_median = median, a_raw_min = min, a_raw_max = max,
                                     models_without_8, a_raw_median_without_8 = median_without_8)], by = "instrument", all.x = TRUE)
F8 <- STAB[instrument == "sf8", .(instrument = "sf8", label = "The 8", stability_models = .N, a_raw_median = median(est), a_raw_min = min(est), a_raw_max = max(est))]
s8 <- rbindlist(lapply(Filter(function(u) !is.null(u) && !is.null(u$a_size) && "sf8" %in% u$a_size$instrument, U_M), function(u) u$a_size[instrument == "sf8"]))
F8[, `:=`(a_std_mean = median(s8$std_mean), a_share_mean = median(s8$share_mean), items = 8L, reversed_items = 3L, response_range = "1-6")]
F_SIZE <- rbind(F8, F_SIZE, fill = TRUE)
fwrite(num(copy(F_SIZE)), sf2b_out("20_f_agreement_size.csv"))
fwrite(num(copy(STAB)), sf2b_out("20_f_agreement_stability.csv"))
# Diagnostic: A free with the measure's factors; beside 04's ri_free correlation of the 8's agreement factor.
RIF04_RI <- fread(sf2b_out("04_agreement_factor_validity.csv"))[form == "sf8", .(score, r_agreement_04_ri_free = r_with_agreement)]
F_DIAG <- rbindlist(lapply(ref_of("a_free_diagnostic"), function(r) {
  if (!r$ok || is.null(r$ra) || !nrow(r$ra)) return(data.table(model = mname(r), error = r$error, improper = r$improper))
  x <- as.data.table(r$ra)[, .(score = sub("^ra_", "", label), r_A_with_factor = est, lower = ci.lower, upper = ci.upper, p = signif(pvalue, 3))]
  rr <- as.data.table(r$est)[grepl("^r_", label) & op == "~~", .(score = sub("^r_", "", label), r_8_in_this_model = est)]
  cbind(model = mname(r), merge(x, rr, by = "score", all.x = TRUE), usable = r$usable, improper = r$improper, error = r$error) }), fill = TRUE)
F_DIAG <- merge(F_DIAG, RIF04_RI, by = "score", all.x = TRUE)
F_DIAG <- merge(F_DIAG, A_OUT[, .(score, r_battery_A)], by = "score", all.x = TRUE)
fwrite(num(copy(F_DIAG)), sf2b_out("20_f_agreement_free_diagnostic.csv"))

## Report ---------
rnd <- function(x) if (is.data.frame(x)) x else if (is.list(x)) lapply(x, rnd) else if (is.numeric(x)) signif(x, 4) else x
gf_json <- function(g) list(tag = g$tag, dropped = g$dropped, smallest_eigenvalue = g$smallest_eigenvalue, smoothed = g$smoothed,
                            parallel_analysis_nfact = g$pa_nfact, group_factors_used = g$K, fa_eigenvalues = g$fa_eigenvalues,
                            g_share_of_common = g$g_share_of_common, g_share_of_total = g$g_share_of_total, second_order_loadings = g$second_order,
                            omega_h = g$omega_h, rms = g$rms, heywood = g$heywood, heywood_without_bfi_intellectual = g$heywood_without_bfi_intellectual,
                            one_factor_share = g$one_factor_share, jennrich_bentler_error = g$jb_error, agreement_between_methods = g$agree,
                            practice_group_net = g$grp_net)
secs_total <- as.numeric(difftime(Sys.time(), T0, units = "secs"))
out <- list(
  analysis = "20: the overlap and the general factor with agreement (acquiescence) modelled across the whole battery. EXPLORATORY, POST HOC, NOT REGISTERED (osf.io/d4a5u does not include it); specified after the confirmatory run and analyses 10-19, with their results seen.",
  mode = SF2B$mode, tag = SF2B$tag, written = format(Sys.time(), "%Y-%m-%d %H:%M"),
  software = list(R = R.version.string, lavaan = as.character(packageVersion("lavaan")), psych = as.character(packageVersion("psych")),
                  GPArotation = as.character(packageVersion("GPArotation"))),
  n = N, n_current = sum(G_DAT$grp == 1), n_never = sum(G_DAT$grp == 0), choices = CHOICES,
  reversed_items = REVTAB, reversed_items_by_score = SCORE_REV,
  checks = c(chk_reg, list(registered_refit_vs_04 = refit_check, ri_free_refit_vs_04 = rif_check, reproduce_13_max_abs_sl_g = repro13)),
  fallbacks = list(table = FALLBACKS, matrix_forms_used = as.list(table(CELLS$form_used)),
                   calibrated_unusable = CELLS[calibrated_source == "primary (calibrated unusable)", model]),
  calibration_loadings = as.list(A_CAL),
  a = list(correlations = A_OUT, h3b = H3B),
  b = list(h3d = H3D),
  c = list(variance_share = C_SHARE, residual_correlations = C_RES),
  d = list(regression = D_REG, residual_correlations = D_RES, practice_group = D_PR),
  e = list(summary = E_SUM, loadings = E_LOAD, net_of_G = E_NET[g_definition == "sl"], practice_group_net = E_GRP,
           details = list(registered_13 = gf_json(GF13), battery_A = gf_json(GFA), calibrated = gf_json(GFC))),
  f = list(size = F_SIZE, free_correlation_diagnostic = F_DIAG),
  models = list(fitted = nrow(FITS), failed = FITS[ok == FALSE, .(part, model, a_mode, error)], improper = FITS[ok == TRUE & nzchar(improper), .(part, model, a_mode, improper)]),
  run_seconds = list(matrix_and_references = round(secs_e1), calibrated = round(secs_cal), parts_b_d = round(secs_cd), bootstrap = round(secs_boot), total = round(secs_total)),
  cores = SF2B$cores)
write_json(rnd(out), sf2b_out("20_exploratory_battery_agreement.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
sf2b_log("20 done in ", round(secs_total), " s; outputs ", sf2b_out("20_*"))
print(num(copy(A_OUT))[, .(score, role, r_registered, r_ri_free, r_battery_A, lower, upper, form_used)])
print(num(copy(H3B))); print(num(copy(H3D))); print(num(copy(C_SHARE))); print(num(copy(C_RES))[, .(score, residual_r_registered, residual_r_battery_A, lower, upper, form_used)])
print(num(copy(D_REG))); print(num(copy(D_RES))[, .(score, residual_r_registered, residual_r_battery_A, lower, upper, form_used)]); print(num(copy(D_PR)))
print(num(copy(E_SUM))); print(num(copy(E_LOAD))[, .(measure, sl_g_13, sl_g_battery_A, sl_g_calibrated, one_factor, sl_g_without_8)])
print(num(copy(E_NET))[g_definition == "sl" & (measure %in% c(FOCUS, EXTENSION)), .(matrix_source, measure, r_with_8, partial_r_net_of_G)]); print(num(copy(E_GRP)))
print(num(copy(F_SIZE))[, .(instrument, form_used, reversed_items, a_raw, a_std_mean, a_share_mean, a_raw_median, a_raw_min, a_raw_max)])
print(num(copy(FALLBACKS))[form_used != "full"])
print(out$checks)
