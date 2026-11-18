# 26: the 8 against FFMQ Observing, criterion by criterion. EXPLORATORY, NOT REGISTERED.
#
# Not part of the filed registration (osf.io/d4a5u) and not in its list of exploratory
# analyses. Specified on 17 November 2026 (placement card S3-22, from the Fable review
# of that date), after the registered analyses and 10-25 had been seen.
#
# Question: on each validity measure and criterion, does the 8 relate differently from
# FFMQ Observing, its nearest neighbour (latent r .865)?
#   latent    one model per measure: 04's model for that measure (the H1 model; the
#             measure's factors, correlated with S and with each other, uncorrelated
#             with the agreement factor) with the FFMQ Observing factor added (for the
#             FFMQ itself it is already there). Labelled: r_<x> = cor(S, x), o_<x> =
#             cor(Observing, x), d_<x> := r_<x> - o_<x>, delta-method 95% CI, Holm
#             over the criteria. Single items (life satisfaction, the past-week items,
#             the engaged pause) as single-indicator factors (04). Brief MAIA-2 with
#             decision 46's bound on maia_12 (as 10, 13, 15).
#   practice  current vs never (H6's groups): a MIMIC with S ~ bS*grp and Observing ~
#             bO*grp (residual variances 1, as 07 and 12), dd := bS - bO.
#   observed  Pearson r of the observed 8 and of observed Observing with each observed
#             score, and the test of two dependent correlations sharing a variable
#             (psych::r.test, Steiger's z with the 8-Observing r).
#
# Settings as registered (Section 5): ordinal items, WLSMV, theta, missing = "pairwise",
# conditional.x = FALSE. Writes aggregates only:
#   output/real/26_latent.csv, 26_practice.csv, 26_observed.csv, 26_fit.csv,
#   26_exploratory_vs_observing.json
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/26_exploratory_vs_observing.R

# Set Up ---------
## Load libraries ---------
packages <- c("jsonlite", "data.table", "lavaan", "psych")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

## Settings: 04's data, nothing fitted (as 12) ---------
source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")
stopifnot(SF2B$mode == "real")
T0 <- Sys.time()
ex04 <- parse("study2b/analysis/04_validity.R", keep.source = FALSE)
first_fit <- which(vapply(ex04, function(e) grepl("^sf2b_log\\(\"fitting", deparse(e)[1]), TRUE))
stopifnot(length(first_fit) == 1)
for (i in seq_len(first_fit - 1)) eval(ex04[[i]], envir = globalenv())
stopifnot(exists("DAT"), exists("score_items"), !exists("RES"))
DAT$grp <- fifelse(M$practice_group == "current", 1L, fifelse(M$practice_group == "never", 0L, NA_integer_))

I8 <- D$sf8; OBS <- "ffmq_observing"; MAIA_BOUND <- 9
SINGLE_CRIT <- c("lifesat", "crit_pause_days", "crit_notice_first", "crit_outdoors", "engaged_pause")
BLOCKS <- SF2B_MEASURES
CHOICES <- list(
  models = "04's model for each measure with the FFMQ Observing factor added (the FFMQ model already holds it). Observing correlates freely with S and with every validity factor, and not with the agreement factor (as every validity factor).",
  maia = "Brief MAIA-2 with decision 46's bound (l_maia12^2 < 9 on maia_12, as 10, 13 and 15).",
  single = "Single items as single-indicator factors with no residual variance (04).",
  practice = "MIMIC on cases with a practice group (current 1, never 0): S ~ bS*grp, Observing ~ bO*grp, residual variances 1; the agreement factor unregressed.",
  holm = "Holm over the latent differences of every criterion except Observing itself (one family).",
  observed = "Observed scores as 01 scored them (A$scores); pairwise complete cases; psych::r.test(n, r12, r13, r23) with r12 = r(8, x), r13 = r(Observing, x), r23 = r(8, Observing).")

## Models ---------
model_of <- function(t) {
  syn <- sf2b_scale_syntax(I8); ord <- I8
  if (!is.null(t$single)) {
    v <- t$single
    syn <- c(syn, sf2b_validity_syntax(OBS, D),
             paste0("L_", v, " =~ 1*", v), paste0(v, " ~~ 0*", v), paste0("L_", v, " ~~ 1*L_", v),
             paste0("S ~~ r_", v, "*L_", v), paste0("RI ~~ 0*L_", v),
             paste0(OBS, " ~~ o_", v, "*L_", v), paste0("d_", v, " := r_", v, " - o_", v))
    ord <- c(ord, score_items[[OBS]], v)
    return(list(syn = syn, ord = unique(ord), crit = v))
  }
  scores <- unique(c(OBS, t$scores)); crit <- setdiff(scores, OBS)
  syn <- c(syn, sf2b_validity_syntax(scores, D))
  if ("maia_notworrying" %in% scores) syn <- c(syn, "maia_notworrying =~ l_maia12*maia_12", sprintf("l_maia12^2 < %s", MAIA_BOUND))
  syn <- c(syn, paste0(OBS, " ~~ o_", crit, "*", crit), paste0("d_", crit, " := r_", crit, " - o_", crit))
  list(syn = syn, ord = unique(c(ord, unlist(score_items[scores]))), crit = crit)
}
TASKS <- c(lapply(names(BLOCKS), function(b) list(block = b, scores = BLOCKS[[b]], single = NULL)),
           lapply(SINGLE_CRIT, function(v) list(block = v, scores = character(), single = v)))
run_task <- function(t) {
  m <- model_of(t)
  f <- sf2b_try(sf2b_sem(m$syn, DAT, ordered = m$ord))
  if (sf2b_failed(f)) return(list(task = t, ok = FALSE, error = sf2b_msg(f)))
  bad <- sf2b_orientation(f, I8, if (is.null(t$single)) unique(c(OBS, t$scores)) else OBS, D)
  if (length(bad)) return(list(task = t, ok = FALSE, error = paste("not oriented as keyed:", paste(bad, collapse = ", "))))
  pe <- as.data.table(parameterEstimates(f))
  out <- rbindlist(lapply(m$crit, function(x) {
    g <- function(l) pe[label == l]
    data.table(block = t$block, criterion = x, r_8 = g(paste0("r_", x))$est, r_8_lo = g(paste0("r_", x))$ci.lower, r_8_hi = g(paste0("r_", x))$ci.upper,
               r_obs = g(paste0("o_", x))$est, r_obs_lo = g(paste0("o_", x))$ci.lower, r_obs_hi = g(paste0("o_", x))$ci.upper,
               diff = g(paste0("d_", x))$est, diff_lo = g(paste0("d_", x))$ci.lower, diff_hi = g(paste0("d_", x))$ci.upper, p = g(paste0("d_", x))$pvalue)
  }))
  r8o <- pe[label == paste0("r_", OBS)]$est
  list(task = t, ok = TRUE, r = out, r8obs = r8o, fit = c(list(block = t$block), sf2b_fit_stats(f)))
}
sf2b_log("fitting ", length(TASKS), " models")
RES <- sf2b_par_lapply(TASKS, run_task, data = list(DAT = DAT, D = D, I8 = I8, OBS = OBS, MAIA_BOUND = MAIA_BOUND,
                                                    score_items = score_items, model_of = model_of))
for (r in Filter(function(r) !r$ok, RES)) sf2b_log("FAILED: ", r$task$block, ": ", r$error)
lat <- rbindlist(lapply(Filter(function(r) r$ok, RES), `[[`, "r"))
lat[, p_holm := p.adjust(p, "holm")]
lat <- merge(lat, SF2B_PREDICTED[, .(criterion = score, label)], by = "criterion", all.x = TRUE, sort = FALSE)
lat[is.na(label), label := criterion]
fits <- rbindlist(lapply(Filter(function(r) r$ok, RES), `[[`, "fit"), fill = TRUE)
r8obs <- vapply(Filter(function(r) r$ok, RES), `[[`, 0, "r8obs")

## Practice MIMIC ---------
pm_syn <- c(sf2b_scale_syntax(I8), sf2b_validity_syntax(OBS, D), "S ~ bS*grp", paste0(OBS, " ~ bO*grp"), "dd := bS - bO")
dp <- DAT[!is.na(DAT$grp), , drop = FALSE]
fp <- sf2b_sem(pm_syn, dp, ordered = unique(c(I8, score_items[[OBS]])))
stopifnot(!sf2b_failed(fp))
pp <- as.data.table(parameterEstimates(fp))[label %in% c("bS", "bO", "dd")]
practice <- pp[, .(quantity = c(bS = "8 (S): current - never", bO = "Observing: current - never", dd = "difference")[label],
                   est, se, ci.lower, ci.upper, pvalue)]
practice[, n := nrow(dp)]

## Observed ---------
SC <- A$scores
obs_rows <- rbindlist(lapply(setdiff(c(unlist(BLOCKS), SINGLE_CRIT), OBS), function(x) {
  ok <- complete.cases(SC$sf8, SC[[OBS]], SC[[x]]); n <- sum(ok)
  r12 <- cor(SC$sf8[ok], SC[[x]][ok]); r13 <- cor(SC[[OBS]][ok], SC[[x]][ok]); r23 <- cor(SC$sf8[ok], SC[[OBS]][ok])
  tt <- psych::r.test(n = n, r12 = r12, r13 = r13, r23 = r23)
  data.table(criterion = x, n = n, r_8 = r12, r_obs = r13, r_8_obs = r23, diff = r12 - r13, z = unname(tt$t), p = unname(tt$p))
}))
obs_rows[, p_holm := p.adjust(p, "holm")]

## Write ---------
num <- function(DT, k = 4) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, function(x) round(x, k)), .SDcols = cols]; DT }
fwrite(num(copy(lat)), sf2b_out("26_latent.csv")); fwrite(num(copy(practice)), sf2b_out("26_practice.csv"))
fwrite(num(copy(obs_rows)), sf2b_out("26_observed.csv")); fwrite(fits, sf2b_out("26_fit.csv"))
write_json(list(
  analysis = "26: the 8 against FFMQ Observing, criterion by criterion. EXPLORATORY, NOT REGISTERED (osf.io/d4a5u does not include it); specified 17 November 2026 after 10-25 were seen (card S3-22).",
  choices = CHOICES, r_8_observing_by_model = round(r8obs, 4), failed = lapply(Filter(function(r) !r$ok, RES), function(r) list(block = r$task$block, error = r$error)),
  latent = num(copy(lat)), practice = num(copy(practice)), observed = num(copy(obs_rows)),
  minutes = round(as.numeric(difftime(Sys.time(), T0, units = "mins")), 1)),
  sf2b_out("26_exploratory_vs_observing.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
print(num(copy(lat))[, .(criterion, r_8, r_obs, diff, diff_lo, diff_hi, p_holm)], nrows = 100)
print(num(copy(practice)))
print(num(copy(obs_rows))[, .(criterion, r_8, r_obs, diff, p_holm)], nrows = 100)
sf2b_log("26 done in ", round(as.numeric(difftime(Sys.time(), T0, units = "mins")), 1), " min")
