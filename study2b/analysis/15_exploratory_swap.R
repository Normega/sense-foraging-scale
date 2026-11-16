# 15: sf6_44 in place of U4. EXPLORATORY, NOT REGISTERED. A POST HOC CHANGE.
#
# Not part of the filed registration (osf.io/d4a5u) and not in its list of exploratory
# analyses. Specified on 15 November 2026, after the confirmatory real-mode run, with
# every loading already seen: U4 (sf2b_07) loads .18 on the 8's content factor, the
# unscored reversed sf6_44 (sf2b_11) .53 in the 13-statement model (03). Choosing the
# replacement from the same data that showed U4's weakness capitalizes on chance, so
# every number here is optimistic for the swapped scale and a change to the scale would
# need a new sample to confirm it. Nothing here changes a registered cell or output, and
# the registered 8 remains the scale.
#
# Three item sets, each run through the same code:
#   registered   the 8 as registered (the reference; should reproduce 03, 04 and 07)
#   swap         the 8 with sf6_44 in place of U4 (sf6_44 reversed, as built)
#   seven        the 8 without U4, and without sf6_44
# For each:
#   H1     the registered model (one content factor + the agreement factor) through the
#          registered rule (sf2b_h1: CFI, the largest residual and the local rule, the
#          close-fit test), RMSEA and SRMR, content omega (03's formula), alpha (03's),
#          the loadings;
#   r      the latent r with the four H2 neighbours (FFMQ Observing from the FFMQ
#          model, SBI, BFI-2 Aesthetic Sensitivity from the BFI-2 model, Brief MAIA-2
#          Noticing from the Brief MAIA-2 model with decision 46's bound on maia_12, the
#          repair that gave the confirmatory cell), NCS-6, both BIDR scales, PSS-4, PHQ-4,
#          life satisfaction, past-week pause days and the engaged pause (single-indicator
#          factors, 04's treatment); every score of each measure's model is written;
#   H6     current vs never on the observed score (07's sf2b_kg: Student's t, pooled-SD
#          d), and its sensitivity rows (spirituality and religiosity; dedicated vs
#          never), and the latent difference (07's MIMIC).
#
# Settings as registered (Section 5): ordinal items, WLSMV, theta, missing = "pairwise",
# conditional.x = FALSE, the agreement factor as registered (loadings 1 on every
# unrecoded statement of the set, variance free, uncorrelated with every other factor).
# CHOICES below are written to the JSON.
#
# Writes output/<tag>/15_*.csv and 15_exploratory_swap.json (aggregates only).
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/15_exploratory_swap.R

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
# 10-14 do): A, D, M, S, DAT, score_items; nothing is fitted.
ex04 <- parse("study2b/analysis/04_validity.R", keep.source = FALSE)
first_fit <- which(vapply(ex04, function(e) grepl("^sf2b_log\\(\"fitting", deparse(e)[1]), TRUE))
stopifnot(length(first_fit) == 1)
for (i in seq_len(first_fit - 1)) eval(ex04[[i]], envir = globalenv())
stopifnot(exists("DAT"), exists("score_items"), !exists("RES"))
ex10 <- parse("study2b/analysis/10_deviation_maia.R", keep.source = FALSE)
is_def <- function(e, nm) is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(nm))
eval(Filter(function(e) is_def(e, "improper_of"), as.list(ex10))[[1]], envir = globalenv())
improper_list <- function(f) { x <- improper_of(f)$list; x[!grepl("^std loading (L_[[:alnum:]_]+|grp) =~ ", x)] }

REV <- D$reversed; SRC <- D$source
stopifnot(SRC[["sf2b_07"]] == "U4", SRC[["sf2b_11"]] == "sf6_44", "sf2b_11" %in% REV, "sf2b_07" %in% REV)
SETS <- list(registered = D$sf8,
             swap = c(setdiff(D$sf8, "sf2b_07"), "sf2b_11"),
             seven = setdiff(D$sf8, "sf2b_07"))
SET_LABEL <- c(registered = "the registered 8", swap = "the 8 with sf6_44 in place of U4", seven = "the 7 (without U4 or sf6_44)")
MAIA_BOUND <- 9
BLOCKS <- list(ffmq = SF2B_MEASURES$ffmq, sbi = "sbi_moment", bfi = SF2B_MEASURES$bfi, maia = SF2B_MEASURES$maia,
               ncs = "ncs6", bidr = SF2B_MEASURES$bidr, pss = "pss4", phq = "phq4_total",
               lifesat = "lifesat", crit_pause_days = "crit_pause_days", engaged_pause = "engaged_pause")
SINGLE1 <- c("lifesat", "crit_pause_days", "engaged_pause")
FOCUS <- c("ffmq_observing", "sbi_moment", "bfi_aesthetic", "maia_noticing", "ncs6", "bidr_sde", "bidr_im", "pss4", "phq4_total",
           "lifesat", "crit_pause_days", "engaged_pause")

CHOICES <- list(
  post_hoc = "sf6_44 was chosen as U4's replacement after its loading in this sample was seen (03's 13-statement model). The swap therefore capitalizes on chance and its numbers are optimistic; confirming it needs a new sample. The registered 8 stays the scale.",
  keying = "sf6_44 is scored reversed, as built (it is reverse-keyed in the scale definition). The agreement factor loads 1 on every statement of each set, unrecoded, as registered.",
  minimum_answered = "Observed scores: 7 of 8 answered for the 8 and the swap (the registered rule), 6 of 7 for the 7 (all but one, the same rule).",
  h1 = "sf2b_h1 (the registered H1 rule and quantities) and 03's alpha (complete cases, reversed statements recoded 7 - x) for each set. Model-specific cutoffs are not computed.",
  validity = "04's model for each measure with the set in place of the 8 (subscales of one measure as correlated factors; validity factors uncorrelated with the agreement factor). Brief MAIA-2 with decision 46's bound (maia_12 standardized loading at most .949, as 10 and 13), so its cells are conditional on the bound; the registered set is refitted the same way, which reproduces 10's Noticing cell. Single items as single-indicator factors with no residual variance (04). Statements are ordered unreversed first, so S takes its sign from an unreversed statement.",
  h6 = "07's sf2b_kg (Student's t, d = coefficient over the pooled SD, d's CI the coefficient's CI over that SD) on each set's observed score: current vs never, adjusted for spirituality and religiosity, dedicated vs never. Latent difference: 07's MIMIC (S ~ group, S's residual variance 1, agreement factor ~ 0*group), cases with a group.",
  improper = "10's improper_of() on every fitted model; single-indicator loadings of 1 by construction are not listed.")

## Models ---------
order_items <- function(items) items[order(items %in% REV)]
block_model <- function(items, b) {
  syn <- sf2b_scale_syntax(items); ord <- items; scores <- character()
  if (b %in% SINGLE1) {
    v <- b
    syn <- c(syn, paste0("L_", v, " =~ 1*", v), paste0(v, " ~~ 0*", v), paste0("L_", v, " ~~ 1*L_", v),
             paste0("S ~~ r_", v, "*L_", v), paste0("RI ~~ 0*L_", v))
    ord <- c(ord, v)
  } else {
    scores <- BLOCKS[[b]]
    syn <- c(syn, sf2b_validity_syntax(scores, D))
    if ("maia_notworrying" %in% scores) syn <- c(syn, "maia_notworrying =~ l_maia12*maia_12", sprintf("l_maia12^2 < %s", MAIA_BOUND))
    ord <- c(ord, unlist(score_items[scores]))
  }
  list(syn = syn, ord = unique(ord), scores = scores)
}
TASKS <- list()
for (s in names(SETS)) {
  TASKS[[length(TASKS) + 1]] <- list(set = s, kind = "h1")
  TASKS[[length(TASKS) + 1]] <- list(set = s, kind = "mimic")
  for (b in names(BLOCKS)) TASKS[[length(TASKS) + 1]] <- list(set = s, kind = "validity", block = b)
}
run_task <- function(t) {
  items <- order_items(SETS[[t$set]]); t0 <- Sys.time()
  base <- list(task = t, ok = FALSE, error = "", improper = "")
  if (t$kind == "h1") {
    h <- sf2b_try(sf2b_h1(X, items, REV))
    if (inherits(h, "error") || !isTRUE(h$converged)) { base$error <- if (inherits(h, "error")) conditionMessage(h) else h$error; return(base) }
    pe <- parameterEstimates(h$fit); riv <- pe[pe$lhs == "RI" & pe$op == "~~" & pe$rhs == "RI", ]
    return(list(task = t, ok = TRUE, error = "", improper = paste(improper_list(h$fit), collapse = "; "),
                h1 = h[setdiff(names(h), c("fit", "residuals"))], residuals = h$residuals))
  }
  if (t$kind == "mimic") {   # 07's mimic(), with the set's statements
    d <- as.data.frame(X[, items]); d$grp <- g_cur; d <- d[!is.na(d$grp), ]
    syn <- c(paste0("S =~ NA*", paste(items, collapse = " + ")), "S ~~ 1*S",
             paste0("RI =~ ", paste0("1*", items, collapse = " + ")), "RI ~~ RI", "S ~~ 0*RI", "S ~ b*grp", "RI ~ 0*grp")
    f <- sf2b_try(sf2b_sem(syn, d, ordered = items))
    if (sf2b_failed(f)) { base$error <- sf2b_msg(f); return(base) }
    pe <- parameterEstimates(f); b <- pe[pe$label == "b", ]
    return(list(task = t, ok = TRUE, error = "", improper = paste(improper_list(f), collapse = "; "),
                mimic = data.table(n = nrow(d), latent_d = b$est, lower = b$ci.lower, upper = b$ci.upper, p = b$pvalue)))
  }
  m <- block_model(items, t$block)
  f <- sf2b_try(sf2b_sem(m$syn, DAT, ordered = m$ord))
  if (sf2b_failed(f)) { base$error <- sf2b_msg(f); return(base) }
  bad <- sf2b_orientation(f, items, m$scores, D)
  pe <- parameterEstimates(f)
  r <- pe[grepl("^r_", pe$label), c("label", "est", "se", "ci.lower", "ci.upper", "pvalue")]
  extra <- list()
  if (t$block == "maia") extra$maia12_std_residual <- diag(lavInspect(f, "std")$theta)[["maia_12"]]
  c(list(task = t, ok = !length(bad), error = if (length(bad)) paste("not oriented as keyed:", paste(bad, collapse = ", ")) else "",
         improper = paste(improper_list(f), collapse = "; "), r = r, fit = sf2b_fit_stats(f),
         secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))), extra)
}
g_cur <- fifelse(M$practice_group == "current", 1L, fifelse(M$practice_group == "never", 0L, NA_integer_))
g_ded <- fifelse(M$dedicated, 1L, fifelse(M$practice_group == "never", 0L, NA_integer_))
sf2b_log("15: ", length(TASKS), " models on ", SF2B$cores, " workers")
t1 <- Sys.time()
RES <- sf2b_par_lapply(TASKS, run_task, data = list(DAT = DAT, X = A$scale_raw, D = D, score_items = score_items, REV = REV, SETS = SETS,
                                                     BLOCKS = BLOCKS, SINGLE1 = SINGLE1, MAIA_BOUND = MAIA_BOUND, g_cur = g_cur,
                                                     order_items = order_items, block_model = block_model,
                                                     improper_of = improper_of, improper_list = improper_list))
secs_models <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
tid <- function(t) paste(t$set, t$kind, t$block %||% "")
for (r in Filter(function(r) !r$ok, RES)) sf2b_log("FAILED: ", tid(r$task), ": ", r$error)
for (r in Filter(function(r) r$ok && nzchar(r$improper), RES)) sf2b_log("IMPROPER: ", tid(r$task), ": ", r$improper)
num <- function(DT, k = 3) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, function(x) round(x, k)), .SDcols = cols]; DT }
get_res <- function(s, kind, b = NULL) Filter(function(r) r$task$set == s && r$task$kind == kind && identical(r$task$block, b), RES)[[1]]

## H1 ---------
alpha_of <- function(items) { Z <- A$scale_raw[, items]; rv <- intersect(items, REV); Z[, rv] <- 7 - Z[, rv]   # 03's
  suppressWarnings(psych::alpha(Z[complete.cases(Z), ], check.keys = FALSE, warnings = FALSE)$total$raw_alpha) }
FIT <- rbindlist(lapply(names(SETS), function(s) {
  r <- get_res(s, "h1")
  if (!r$ok) return(data.table(set = s, ok = FALSE, error = r$error))
  h <- r$h1; st <- h$stats
  data.table(set = s, label = SET_LABEL[[s]], statements = paste(SRC[order_items(SETS[[s]])], collapse = " "), ok = TRUE, error = "",
             n = st$ntotal, chisq = st$chisq.scaled, df = st$df, p = st$pvalue.scaled, cfi = st$cfi.scaled, tli = st$tli.scaled,
             rmsea = st$rmsea.scaled, rmsea_lower = st$rmsea.ci.lower.scaled, rmsea_upper = st$rmsea.ci.upper.scaled, srmr = st$srmr,
             max_residual = h$max_residual, max_residual_pair = h$max_residual_pair,
             max_residual_pair_ids = paste(SRC[strsplit(h$max_residual_pair, " ~~ ")[[1]]], collapse = " ~~ "),
             freeing_it_changes_omega_by = h$d_omega, freeing_it_changes_a_loading_by_up_to = h$d_lambda,
             usrmr = h$usrmr, closefit_p = h$closefit_p, cfi_ok = h$cfi_ok, local_ok = h$local_ok, close_ok = h$close_ok,
             h1_rule_met = h$supported, omega = h$omega, alpha = alpha_of(SETS[[s]]),
             agreement_variance = h$ri_variance, agreement_share = h$ri_share, improper = r$improper)
}), fill = TRUE)
fwrite(num(copy(FIT), 4), sf2b_out("15_fit.csv"))
LOAD <- rbindlist(lapply(names(SETS), function(s) { r <- get_res(s, "h1"); if (!r$ok) return(NULL)
  cbind(set = s, r$h1$loadings[, .(statement = item, source = SRC[item], reversed = item %in% REV, content, keyed_content, agreement, residual_variance)]) }))
fwrite(num(copy(LOAD)), sf2b_out("15_loadings.csv"))
RESID <- rbindlist(lapply(names(SETS), function(s) { r <- get_res(s, "h1"); if (!r$ok) return(NULL)
  R <- r$residuals; it <- rownames(R); ij <- which(upper.tri(R), arr.ind = TRUE)
  data.table(set = s, a = SRC[it[ij[, 1]]], b = SRC[it[ij[, 2]]], residual = R[ij]) }))
RESID <- RESID[order(set, -abs(residual))]
fwrite(num(copy(RESID)), sf2b_out("15_residuals.csv"))

## Validity ---------
VAL <- rbindlist(lapply(Filter(function(r) r$task$kind == "validity", RES), function(r) {
  if (!r$ok) return(data.table(set = r$task$set, block = r$task$block, ok = FALSE, error = r$error))
  data.table(set = r$task$set, block = r$task$block, ok = TRUE, error = "", score = sub("^r_", "", r$r$label), r = r$r$est, se = r$r$se,
             lower = r$r$ci.lower, upper = r$r$ci.upper, p = r$r$pvalue, improper = r$improper,
             maia12_std_residual = r$maia12_std_residual %||% NA_real_, cfi = r$fit$cfi.scaled, rmsea = r$fit$rmsea.scaled, srmr = r$fit$srmr)
}), fill = TRUE)
VAL <- merge(VAL, SF2B_PREDICTED[, .(score, label)], by = "score", all.x = TRUE, sort = FALSE)
VAL[, focus := score %in% FOCUS]
VAL <- VAL[order(match(set, names(SETS)), -focus, match(score, SF2B_PREDICTED$score))]
fwrite(num(copy(VAL)), sf2b_out("15_validity.csv"))
# The registered set against 04 (and 10 for Noticing).
LAT04 <- fread(sf2b_out("04_latent_validity.csv"))[form == "sf8" & sample == "all" & ri_free == FALSE]
chk <- merge(VAL[set == "registered" & ok == TRUE, .(score, r_here = r)], LAT04[, .(score, r_04 = est)], by = "score", all.x = TRUE)
n10 <- tryCatch(fromJSON(sf2b_out("10_deviation_maia.json"))$sf8$h2_noticing$r, error = function(e) NA_real_)
chk[score == "maia_noticing", r_04 := n10]
chk[, diff := round(r_here - r_04, 3)]

## H6 ---------
SCX <- A$scale_raw
KG <- rbindlist(lapply(names(SETS), function(s) {
  it <- SETS[[s]]; y <- sf2b_score(SCX, it, REV, need = length(it) - 1L)
  rbind(cbind(set = s, analysis = "H6: current vs never", sf2b_kg(y, g_cur)),
        cbind(set = s, analysis = "adjusted for spirituality and religiosity", sf2b_kg(y, g_cur, M[, .(spiritual, religious)])),
        cbind(set = s, analysis = "dedicated vs never", sf2b_kg(y, g_ded)))
}))
MIM <- rbindlist(lapply(names(SETS), function(s) { r <- get_res(s, "mimic")
  if (!r$ok) data.table(set = s, error = r$error) else cbind(set = s, r$mimic, improper = r$improper) }), fill = TRUE)
OBS <- data.table(set = names(SETS), n_scored = vapply(SETS, function(it) sum(!is.na(sf2b_score(SCX, it, REV, need = length(it) - 1L))), 0),
                  r_with_registered_8 = vapply(SETS, function(it) cor(sf2b_score(SCX, it, REV, need = length(it) - 1L), S$sf8, use = "complete.obs"), 0))
fwrite(num(copy(KG), 4), sf2b_out("15_known_groups.csv"))
fwrite(num(copy(MIM), 4), sf2b_out("15_latent_difference.csv"))
chk07 <- fread(sf2b_out("07_known_groups.csv"))[form == "sf8", .(analysis, d_07 = d)]
chk07 <- merge(chk07, KG[set == "registered", .(analysis, d_here = d)], by = "analysis")

## Side by side ---------
wide <- function(DT, val, by) { x <- dcast(DT, as.formula(paste(by, "~ set")), value.var = val); setcolorder(x, c(by, names(SETS))); x }
CMP <- rbind(
  data.table(quantity = c("CFI", "RMSEA", "SRMR", "largest |residual|", "H1 rule met", "content omega", "alpha", "agreement share"),
             registered = NA_character_, swap = NA_character_, seven = NA_character_),
  fill = TRUE)
fmt <- function(x) if (is.logical(x)) as.character(x) else sprintf("%.3f", x)
for (s in names(SETS)) { f <- FIT[set == s]
  CMP[[s]] <- c(fmt(f$cfi), fmt(f$rmsea), fmt(f$srmr), fmt(f$max_residual), fmt(f$h1_rule_met), fmt(f$omega), fmt(f$alpha), fmt(f$agreement_share)) }
vr <- VAL[focus == TRUE & ok == TRUE, .(quantity = paste("latent r:", score), set, v = sprintf("%.3f [%.3f, %.3f]", r, lower, upper))]
vr <- dcast(vr, quantity ~ set, value.var = "v")[order(match(sub("latent r: ", "", quantity), FOCUS))]
kr <- KG[, .(quantity = paste("d:", analysis), set, v = sprintf("%.3f [%.3f, %.3f]", d, d_lower, d_upper))]
kr <- dcast(kr, quantity ~ set, value.var = "v")
mr <- MIM[, .(quantity = "latent d (MIMIC)", set, v = sprintf("%.3f [%.3f, %.3f]", latent_d, lower, upper))]
mr <- dcast(mr, quantity ~ set, value.var = "v")
CMP <- rbind(CMP, vr, kr, mr, fill = TRUE)
setcolorder(CMP, c("quantity", names(SETS)))
fwrite(CMP, sf2b_out("15_side_by_side.csv"))

## Report ---------
rnd <- function(x) if (is.data.frame(x)) x else if (is.list(x)) lapply(x, rnd) else if (is.numeric(x)) signif(x, 4) else x
out <- list(
  analysis = "15: sf6_44 in place of U4. EXPLORATORY, NOT REGISTERED (osf.io/d4a5u does not include it). A POST HOC CHANGE: the replacement was chosen after this sample's loadings were seen, so its numbers are optimistic and any change to the scale needs a new sample to confirm it. The registered 8 remains the scale.",
  mode = SF2B$mode, tag = SF2B$tag, written = format(Sys.time(), "%Y-%m-%d %H:%M"),
  software = list(R = R.version.string, lavaan = as.character(packageVersion("lavaan")), psych = as.character(packageVersion("psych"))),
  n = nrow(DAT), choices = CHOICES, sets = lapply(SETS, function(it) unname(SRC[order_items(it)])),
  models = list(fitted = length(RES), failed = lapply(Filter(function(r) !r$ok, RES), function(r) list(model = tid(r$task), error = r$error)),
                improper = lapply(Filter(function(r) r$ok && nzchar(r$improper), RES), function(r) list(model = tid(r$task), improper = r$improper)),
                maia12_std_residual = VAL[block == "maia" & ok == TRUE, .(maia12_std_residual = maia12_std_residual[1]), by = set]),
  reference_checks = list(validity_against_04_and_10 = chk, known_groups_against_07 = chk07),
  side_by_side = CMP, fit = FIT, loadings = LOAD, validity_focus = VAL[focus == TRUE], known_groups = KG, latent_difference = MIM,
  observed_scores = OBS,
  run_seconds = list(models = round(secs_models), total = round(as.numeric(difftime(Sys.time(), T0, units = "secs")))), cores = SF2B$cores)
write_json(rnd(out), sf2b_out("15_exploratory_swap.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
sf2b_log("15 done in ", out$run_seconds$total, " s; outputs ", sf2b_out("15_*"))
print(chk); print(chk07)
print(CMP)
print(num(copy(LOAD))[, .(set, source, keyed_content, agreement)])
print(OBS)
