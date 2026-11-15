# 12: what is left of the 8 once FFMQ Observing is accounted for. EXPLORATORY, NOT REGISTERED.
#
# Not part of the filed registration (osf.io/d4a5u) and not in its list of exploratory
# analyses. Run after the confirmatory real-mode run, with its answers seen; nothing here
# changes a registered cell or output. Every number it writes is exploratory.
#
# Why: the 8's latent correlation with FFMQ Observing is .865 (04, the FFMQ model), the
# H3a upper bound .912. This asks what the 8 measures that Observing does not.
#   (a) Residual model. The H1 model plus the FFMQ Observing factor, with the 8's content
#       factor S regressed on Observing (S ~ b*Obs, S's residual variance 1). S's
#       disturbance is the part of the 8 not shared with Observing. One model per
#       validity measure adds that measure's factor(s), correlated with Observing and with
#       S's disturbance: the disturbance's correlations with the measures. The mirror
#       model regresses Observing on S, for Observing's own residual. Practice group
#       (current vs never) by MIMIC; past-week pause days as a single-indicator factor.
#   (b) Bifactor check over the 8's and Observing's items, with the agreement factor.
#   (c) Item level: cross-measure polychorics and residual correlations of the two-factor
#       model. Item ids only: FFMQ wording is copyrighted and is never written.
#
# Settings as registered (Section 5): ordinal items, WLSMV, theta, missing = "pairwise",
# conditional.x = FALSE, the 8's agreement factor (loadings 1 on the unrecoded statements,
# variance free, uncorrelated with every other factor). The choices this script makes
# where the registration says nothing are listed in CHOICES below and written to the JSON.
#
# Writes output/<tag>/12_*.csv and 12_exploratory_overlap.json (aggregates only).
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/12_exploratory_overlap.R

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
T0 <- Sys.time()

# The pipeline's own data: 04_validity.R's top-level expressions up to (not including)
# its first fit (as 10_deviation_maia.R does). They build A, D, M, S, DAT, score_items
# and fit nothing.
ex04 <- parse("study2b/analysis/04_validity.R", keep.source = FALSE)
first_fit <- which(vapply(ex04, function(e) grepl("^sf2b_log\\(\"fitting", deparse(e)[1]), TRUE))
stopifnot(length(first_fit) == 1)
for (i in seq_len(first_fit - 1)) eval(ex04[[i]], envir = globalenv())
stopifnot(exists("DAT"), exists("score_items"), !exists("RES"))
# 10_deviation_maia.R's check for improper solutions, taken as written.
ex10 <- parse("study2b/analysis/10_deviation_maia.R", keep.source = FALSE)
is_def <- function(e, nm) is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(nm))
eval(Filter(function(e) is_def(e, "improper_of"), as.list(ex10))[[1]], envir = globalenv())
# A single-indicator factor (no residual variance) and lavaan's phantom for an observed
# covariate have a standardized loading of 1 by construction: not improper, so not listed.
improper_list <- function(f) { x <- improper_of(f)$list; x[!grepl("^std loading (L_[[:alnum:]_]+|grp) =~ ", x)] }

# Practice group as 07_known_groups.R codes it: current = 1, never = 0, everyone else NA.
DAT$grp <- fifelse(M$practice_group == "current", 1L, fifelse(M$practice_group == "never", 0L, NA_integer_))

OBS <- "ffmq_observing"
I8 <- D$sf8
IO <- score_items[[OBS]]
REV8 <- intersect(I8, D$reversed)

# Targets: the other three H2 neighbours, NCS-6, both BIDR scales (one model, as their
# measure), PSS-4, PHQ-4, life satisfaction, past-week pause days and the engaged pause.
# Added (not asked for, cheap, and named in 13): EQ Decentering, SAMD, BFI-2 Creative.
TARGETS <- list(sbi_moment = "sbi_moment", bfi_aesthetic = "bfi_aesthetic", maia_noticing = "maia_noticing",
                ncs6 = "ncs6", bidr = c("bidr_sde", "bidr_im"), pss4 = "pss4", phq4_total = "phq4_total",
                lifesat = "lifesat", crit_pause_days = "crit_pause_days", engaged_pause = "engaged_pause",
                eq_decentering = "eq_decentering", samd_md = "samd_md", bfi_creative = "bfi_creative")
SINGLE <- c("lifesat", "crit_pause_days", "engaged_pause")

CHOICES <- list(
  residual_definition = "S's disturbance in S ~ b*Obs (S's residual variance fixed at 1; Obs variance 1), so 1 - R^2 = 1/(1 + b^2). The mirror model regresses Obs on S the same way.",
  one_model_per_measure = "Each validity measure is added to the residual model in a model of its own (as 04 does), correlated with Obs and with the focal factor's disturbance; the structural part is saturated, so each model's fit equals that of the model with S, Obs and the measure simply correlated. The disturbance's correlation with a measure is therefore the semipartial correlation (measure not residualized).",
  observing_factor = "FFMQ Observing as one factor of its three items (ffmq_01, ffmq_06, ffmq_11), not inside the five-facet FFMQ model; its latent r with the 8 here is reported beside 04's .865.",
  maia_noticing = "Brief MAIA-2 Noticing as one factor of its three items (maia_02-04), without the other subscales, so maia_12's improper loading (decision 46) does not enter.",
  bidr = "BIDR SDE and IM together in one model, as correlated factors (their registered treatment).",
  single_items = "Life satisfaction, past-week pause days and the engaged pause as single-indicator factors with no residual variance, as 04.",
  agreement_factor = "On the 8 only, as registered; uncorrelated with Obs, with every validity factor and (MIMIC) with practice group (RI ~ 0*grp).",
  practice_group = "Current (1) vs never (0) practitioners as 07 codes them; MIMIC as 07 (observed exogenous covariate, conditional.x = FALSE), cases with a group only. S ~ b*Obs + g*grp, Obs ~ a*grp: g is the group difference in the 8's residual in within-group SDs of that residual (d_residual); g/sqrt(1 + b^2) expresses it in within-group SDs of the 8's whole content factor.",
  bifactor = "General factor over the 11 items, a specific factor for the 8, and one for Observing when identified; all orthogonal; agreement factor on the 8 as registered. Loadings start at +/-.5 by keying. Omega counts the agreement factor as error (the pipeline's content omega). If the two-specific model fails or is improper, the bifactor-(S-1) model (Observing as the reference facet, no specific factor for it) is reported.",
  item_level = "Polychorics: lavCor, ordered, pairwise. Keyed sign (the 8's reversed statements sign-flipped) is used for ranking. Residuals: type = cor.bentler (observed minus model-implied polychoric), as H1, in the model of S + RI + Obs correlated.",
  vcov_warning = "Every MIMIC (here and in 07's registered one, checked on the same data) draws lavaan's warning that the parameter covariance matrix has a smallest eigenvalue of about 1e-18: machine zero, from the observed covariate's own parameters under conditional.x = FALSE. It is recorded in the outputs and not treated as a failure.",
  observed_analogues = "Observed-score analogues (not asked for): the observed 8 residualized on observed Observing (OLS) against pause days (Spearman), and the 8 adjusted for observed Observing in the known-groups regression (sf2b_kg).")

# A copy of sf2b_sem that does not swallow warnings, so identification and
# non-positive-definite warnings can be recorded (the call is otherwise identical).
sem_txt <- deparse(sf2b_sem)
stopifnot(sum(grepl("suppressWarnings(do.call(lavaan::sem, args))", sem_txt, fixed = TRUE)) == 1)
sf2b_sem_loud <- eval(parse(text = sub("suppressWarnings(do.call(lavaan::sem, args))", "do.call(lavaan::sem, args)", sem_txt, fixed = TRUE)))
environment(sf2b_sem_loud) <- environment(sf2b_sem)
fit_w <- function(syn, d, ord) {
  w <- character()
  f <- withCallingHandlers(sf2b_try(sf2b_sem_loud(syn, d, ordered = ord)),
                           warning = function(m) { w <<- c(w, conditionMessage(m)); invokeRestart("muffleWarning") })
  list(fit = f, warnings = unique(trimws(gsub("[[:space:]]+", " ", w))))
}

## Syntax helpers ---------
# A validity score's factor (pipeline syntax) without its S ~~ line, keeping RI ~~ 0*.
fac_lines <- function(sc) { x <- sf2b_validity_syntax(sc, D); x[!grepl("^S ~~", x)] }
single_lines <- function(v) c(paste0("L_", v, " =~ 1*", v), paste0(v, " ~~ 0*", v), paste0("L_", v, " ~~ 1*L_", v), paste0("RI ~~ 0*L_", v))
target_factors <- function(tg) ifelse(tg %in% SINGLE, paste0("L_", tg), tg)
# direction "8|obs": focal S regressed on Obs; "obs|8": focal Obs regressed on S.
residual_model <- function(tg, direction) {
  focal <- if (direction == "8|obs") "S" else OBS; other <- if (direction == "8|obs") OBS else "S"
  syn <- c(sf2b_scale_syntax(I8), fac_lines(OBS), paste0(focal, " ~ b*", other),
           "r2 := b^2 / (1 + b^2)", "unexplained := 1 / (1 + b^2)", "r_8_obs := b / sqrt(1 + b^2)")
  ord <- c(I8, IO)
  for (t in tg) {
    fn <- target_factors(t)
    syn <- c(syn, if (t %in% SINGLE) single_lines(t) else fac_lines(t),
             paste0(focal, " ~~ res_", t, "*", fn), paste0(other, " ~~ oth_", t, "*", fn),
             sprintf("tot_%s := (b*oth_%s + res_%s) / sqrt(1 + b^2)", t, t, t))
    ord <- c(ord, if (t %in% SINGLE) t else score_items[[t]])
  }
  list(syn = syn, ord = unique(ord))
}

## (a) Residual models ---------
TASKS <- c(list(list(name = "base", tg = character(), direction = "8|obs"), list(name = "base", tg = character(), direction = "obs|8")),
           unlist(lapply(names(TARGETS), function(n) lapply(c("8|obs", "obs|8"), function(dr) list(name = n, tg = TARGETS[[n]], direction = dr))), recursive = FALSE))
run_res <- function(t) {
  m <- residual_model(t$tg, t$direction)
  t0 <- Sys.time(); fw <- fit_w(m$syn, DAT, m$ord); f <- fw$fit
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (sf2b_failed(f)) return(list(task = t, ok = FALSE, error = sf2b_msg(f), warnings = fw$warnings, secs = secs))
  bad <- sf2b_orientation(f, I8, c(OBS, setdiff(t$tg, SINGLE)), D)
  pe <- parameterEstimates(f)
  keep <- pe$label %in% c("b", "r2", "unexplained", "r_8_obs") | grepl("^(res|oth|tot)_", pe$label)
  list(task = t, ok = !length(bad), error = if (length(bad)) paste("not oriented as keyed:", paste(bad, collapse = ", ")) else "",
       warnings = fw$warnings, secs = secs, est = pe[keep, c("label", "est", "se", "ci.lower", "ci.upper", "pvalue")],
       fit = sf2b_fit_stats(f), improper = paste(improper_list(f), collapse = "; "))
}
sf2b_log("12(a): ", length(TASKS), " residual models on ", SF2B$cores, " workers")
t_a <- Sys.time()
RES <- sf2b_par_lapply(TASKS, run_res, data = list(DAT = DAT, D = D, score_items = score_items, I8 = I8, IO = IO, OBS = OBS,
                                                   SINGLE = SINGLE, residual_model = residual_model, fac_lines = fac_lines,
                                                   single_lines = single_lines, target_factors = target_factors, fit_w = fit_w,
                                                   sf2b_sem_loud = sf2b_sem_loud, improper_of = improper_of, improper_list = improper_list))
secs_a <- as.numeric(difftime(Sys.time(), t_a, units = "secs"))
for (r in Filter(function(r) !r$ok, RES)) sf2b_log("FAILED: ", r$task$name, " ", r$task$direction, ": ", r$error)

fits_a <- rbindlist(lapply(RES, function(r) c(list(part = "a", model = paste(r$task$name, r$task$direction), ok = r$ok, error = r$error,
                                                    improper = r$improper %||% "", warnings = paste(r$warnings, collapse = " | "),
                                                    seconds = round(r$secs, 1)), if (r$ok) r$fit)), fill = TRUE)
EST <- rbindlist(lapply(Filter(function(r) r$ok, RES), function(r) cbind(data.table(target = r$task$name, direction = r$task$direction), as.data.table(r$est))))
EST[, score := sub("^(res|oth|tot)_", "", label)]

base <- EST[target == "base"]
share <- base[label %in% c("b", "r2", "unexplained", "r_8_obs"), .(direction, label, est = r3(est), se = r3(se), lower = r3(ci.lower), upper = r3(ci.upper))]
fwrite(share, sf2b_out("12_a_variance_share.csv"))

labs <- SF2B_PREDICTED[, .(score, label)]
resid <- EST[grepl("^(res|tot)_", label) & target != "base"]
resid[, what := fifelse(grepl("^res_", label), "residual", "total")]
resid <- dcast(resid, score + direction ~ what, value.var = c("est", "ci.lower", "ci.upper", "pvalue"))
resid <- merge(resid, labs, by = "score", all.x = TRUE, sort = FALSE)
resid[, focal := fifelse(direction == "8|obs", "the 8, residual after FFMQ Observing", "FFMQ Observing, residual after the 8")]
LAT04 <- fread(sf2b_out("04_latent_validity.csv"))[form == "sf8" & sample == "all" & ri_free == FALSE]
reg04 <- LAT04[, .(score, r_04 = est)]
resid <- merge(resid, reg04, by = "score", all.x = TRUE, sort = FALSE)
num <- function(DT, k = 3) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, function(x) round(x, k)), .SDcols = cols]; DT }
resid_out <- num(resid[, .(focal, score, label, r_total = est_total, residual_r = est_residual, lower = ci.lower_residual, upper = ci.upper_residual,
                           p = signif(pvalue_residual, 3), r_8_in_04 = r_04)])
resid_out <- resid_out[order(focal, match(score, unlist(TARGETS)))]
fwrite(resid_out, sf2b_out("12_a_residual_correlations.csv"))

## (a) Practice group: MIMIC ---------
G <- DAT[!is.na(DAT$grp), , drop = FALSE]
mimic_syn <- function(kind) {
  core <- c(sf2b_scale_syntax(I8), fac_lines(OBS), "RI ~ 0*grp")
  switch(kind,
         raw = c(core, "S ~ d_8*grp", paste0(OBS, " ~ d_obs*grp"), paste0("S ~~ c_res*", OBS)),
         `8|obs` = c(core, paste0("S ~ b*", OBS, " + g*grp"), paste0(OBS, " ~ a*grp"),
                     "d_residual := g", "d_in_8_sd := g / sqrt(1 + b^2)"),
         `obs|8` = c(core, paste0(OBS, " ~ b*S + g*grp"), "S ~ a*grp",
                     "d_residual := g", "d_in_obs_sd := g / sqrt(1 + b^2)"))
}
mim <- rbindlist(lapply(c("raw", "8|obs", "obs|8"), function(k) {
  fw <- fit_w(mimic_syn(k), G, c(I8, IO)); f <- fw$fit
  if (sf2b_failed(f)) return(data.table(model = k, error = sf2b_msg(f)))
  pe <- parameterEstimates(f); pe <- pe[pe$label %in% c("d_8", "d_obs", "g", "a", "b", "d_residual", "d_in_8_sd", "d_in_obs_sd"), ]
  data.table(model = k, n = nrow(G), n_current = sum(G$grp == 1), n_never = sum(G$grp == 0), label = pe$label, est = pe$est, se = pe$se,
             lower = pe$ci.lower, upper = pe$ci.upper, p = pe$pvalue, error = "", improper = paste(improper_list(f), collapse = "; "),
             warnings = paste(fw$warnings, collapse = " | "))
}), fill = TRUE)
# Observed-score analogues: sf2b_kg (07's function), the 8 adjusted for observed Observing, and the reverse.
g_cur <- DAT$grp
kg_obs <- rbind(cbind(analysis = "observed 8, unadjusted", sf2b_kg(S$sf8, g_cur)),
                cbind(analysis = "observed 8, adjusted for observed FFMQ Observing", sf2b_kg(S$sf8, g_cur, data.table(obs = S$ffmq_observing))),
                cbind(analysis = "observed FFMQ Observing, unadjusted", sf2b_kg(S$ffmq_observing, g_cur)),
                cbind(analysis = "observed FFMQ Observing, adjusted for the observed 8", sf2b_kg(S$ffmq_observing, g_cur, data.table(sf8 = S$sf8))))
fwrite(num(copy(mim)), sf2b_out("12_a_known_groups_latent.csv"))
fwrite(num(copy(kg_obs)), sf2b_out("12_a_known_groups_observed.csv"))

## (a) Past-week pause days: observed analogue ---------
ok <- complete.cases(S$sf8, S$ffmq_observing)
res8 <- rep(NA_real_, nrow(S)); res8[ok] <- resid(lm(sf8 ~ ffmq_observing, data = S[ok]))
resO <- rep(NA_real_, nrow(S)); resO[ok] <- resid(lm(ffmq_observing ~ sf8, data = S[ok]))
pause_obs <- rbindlist(lapply(list(`observed 8` = S$sf8, `observed 8, residual after observed Observing` = res8,
                                   `observed Observing` = S$ffmq_observing, `observed Observing, residual after the observed 8` = resO),
                              function(x) as.data.table(sf2b_r_ci(x, S$crit_pause_days, "spearman"))), idcol = "score")
fwrite(num(copy(pause_obs)), sf2b_out("12_a_pause_observed.csv"))

## (b) Bifactor ---------
keyed_start <- function(it, v = .5) paste0("start(", ifelse(it %in% REV8, -v, v), ")*", it, collapse = " + ")
bif_syn <- function(spec_obs, extra = character()) {
  all_items <- c(I8, IO)
  s <- c(paste0("G =~ NA*", all_items[1], " + ", keyed_start(all_items)), "G ~~ 1*G",
         paste0("S8 =~ NA*", I8[1], " + ", keyed_start(I8)), "S8 ~~ 1*S8",
         paste0("RI =~ ", paste0("1*", I8, collapse = " + ")), "RI ~~ RI",
         "G ~~ 0*S8", "G ~~ 0*RI", "S8 ~~ 0*RI")
  if (spec_obs) s <- c(s, paste0("SO =~ NA*", IO[1], " + ", keyed_start(IO)), "SO ~~ 1*SO", "SO ~~ 0*G", "SO ~~ 0*S8", "SO ~~ 0*RI")
  c(s, extra)
}
bif_indices <- function(f, spec_obs) {
  std <- lavInspect(f, "std"); L <- std$lambda; Th <- std$theta
  all_items <- c(I8, IO); k <- ifelse(all_items %in% REV8, -1, 1); names(k) <- all_items
  col <- function(fac, it) if (fac %in% colnames(L)) L[it, fac] * k[it] else rep(0, length(it))
  th_sum <- function(it) sum(Th[it, it] * outer(k[it], k[it]))
  V_all <- sum(col("G", all_items))^2 + sum(col("S8", all_items))^2 + sum(col("SO", all_items))^2 + sum(col("RI", all_items))^2 + th_sum(all_items)
  V8 <- sum(col("G", I8))^2 + sum(col("S8", I8))^2 + sum(col("RI", I8))^2 + th_sum(I8)
  VO <- sum(col("G", IO))^2 + sum(col("SO", IO))^2 + th_sum(IO)
  list(omega_total_all = (sum(col("G", all_items))^2 + sum(col("S8", all_items))^2 + sum(col("SO", all_items))^2) / V_all,
       omega_h_general_all = sum(col("G", all_items))^2 / V_all,
       ecv_general_all = sum(col("G", all_items)^2) / sum(col("G", all_items)^2 + col("S8", all_items)^2 + col("SO", all_items)^2),
       omega_8 = (sum(col("G", I8))^2 + sum(col("S8", I8))^2) / V8,
       omega_h_general_in_8 = sum(col("G", I8))^2 / V8,
       omega_s_specific_8 = sum(col("S8", I8))^2 / V8,
       specific_share_of_reliable_8 = sum(col("S8", I8))^2 / (sum(col("G", I8))^2 + sum(col("S8", I8))^2),
       ecv_specific_8 = sum(col("S8", I8)^2) / sum(col("G", I8)^2 + col("S8", I8)^2),
       agreement_share_8 = sum(col("RI", I8))^2 / V8,
       omega_obs = (sum(col("G", IO))^2 + sum(col("SO", IO))^2) / VO,
       omega_s_specific_obs = if (spec_obs) sum(col("SO", IO))^2 / VO else NA_real_)
}
bif_run <- function(spec_obs) {
  label <- if (spec_obs) "bifactor: G + S8 + SO" else "bifactor-(S-1): G + S8, Observing as reference"
  t0 <- Sys.time(); fw <- fit_w(bif_syn(spec_obs), DAT, c(I8, IO)); f <- fw$fit
  out <- list(model = label, converged = !sf2b_failed(f), error = sf2b_msg(f), warnings = fw$warnings,
              seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
  if (!out$converged) return(out)
  imp <- improper_of(f)
  vc_ok <- tryCatch(min(eigen(lavInspect(f, "vcov"), symmetric = TRUE, only.values = TRUE)$values) > 0, error = function(e) FALSE)
  std <- lavInspect(f, "std")$lambda
  out <- c(out, list(fit = sf2b_fit_stats(f), improper = imp$list, post_check = imp$post_check, vcov_positive_definite = vc_ok,
                     indices = bif_indices(f, spec_obs),
                     loadings = data.table(item = rownames(std), as.data.table(std))))
  out$proper <- out$converged && !length(imp$list) && vc_ok && !any(grepl("not be identified|not positive definite|could not", fw$warnings))
  # Relations of the specific factor: practice group (MIMIC, cases with a group) and pause days.
  facs <- c("G", "S8", if (spec_obs) "SO")
  mim_extra <- c(paste0(facs, " ~ d_", facs, "*grp"), "RI ~ 0*grp")
  fm <- fit_w(bif_syn(spec_obs, mim_extra), G, c(I8, IO))
  pause_extra <- c(single_lines("crit_pause_days"), paste0(facs, " ~~ p_", facs, "*L_crit_pause_days"))
  fp <- fit_w(bif_syn(spec_obs, pause_extra), DAT, c(I8, IO, "crit_pause_days"))
  rel <- list()
  for (x in list(list(f = fm, what = "practice group (MIMIC d; factor residual variance 1)", lab = "^d_"),
                 list(f = fp, what = "past-week pause days (latent r)", lab = "^p_"))) {
    if (sf2b_failed(x$f$fit)) { rel[[length(rel) + 1]] <- data.table(relation = x$what, error = sf2b_msg(x$f$fit)); next }
    pe <- parameterEstimates(x$f$fit); pe <- pe[grepl(x$lab, pe$label), ]
    rel[[length(rel) + 1]] <- data.table(relation = x$what, factor = sub(x$lab, "", pe$label), est = pe$est, se = pe$se,
                                         lower = pe$ci.lower, upper = pe$ci.upper, p = pe$pvalue, error = "",
                                         improper = paste(improper_list(x$f$fit), collapse = "; "),
                                         warnings = paste(x$f$warnings, collapse = " | "))
  }
  out$relations <- rbindlist(rel, fill = TRUE)
  out
}
sf2b_log("12(b): bifactor models")
t_b <- Sys.time()
BIF <- list(full = bif_run(TRUE), s1 = bif_run(FALSE))
secs_b <- as.numeric(difftime(Sys.time(), t_b, units = "secs"))
for (b in BIF) sf2b_log("  ", b$model, ": converged ", b$converged, "; proper ", isTRUE(b$proper),
                        if (length(b$improper)) paste0(" (", paste(b$improper, collapse = "; "), ")"))
bif_fit <- rbindlist(lapply(BIF, function(b) c(list(model = b$model, converged = b$converged, proper = isTRUE(b$proper),
                                                    vcov_positive_definite = b$vcov_positive_definite %||% NA,
                                                    improper = paste(b$improper %||% "", collapse = "; "), warnings = paste(b$warnings, collapse = " | "),
                                                    error = b$error, seconds = b$seconds), b$fit, b$indices)), fill = TRUE)
fwrite(num(copy(bif_fit), 4), sf2b_out("12_b_bifactor_fit.csv"))
fwrite(num(rbindlist(lapply(BIF, function(b) if (!is.null(b$loadings)) cbind(model = b$model, b$loadings)), fill = TRUE)), sf2b_out("12_b_bifactor_loadings.csv"))
fwrite(num(rbindlist(lapply(BIF, function(b) if (!is.null(b$relations)) cbind(model = b$model, b$relations)), fill = TRUE)), sf2b_out("12_b_bifactor_relations.csv"))

## (c) Item level ---------
sf2b_log("12(c): item level")
t_c <- Sys.time()
PC <- lavCor(as.data.frame(DAT[, c(I8, IO)]), ordered = c(I8, IO), missing = "pairwise")
kk <- ifelse(c(I8, IO) %in% REV8, -1, 1); names(kk) <- c(I8, IO)
cross <- CJ(item_8 = I8, item_obs = IO)
cross[, polychoric := PC[cbind(item_8, item_obs)]][, keyed_polychoric := polychoric * kk[item_8]]
cross[, reversed_in_8 := item_8 %in% REV8]
within8 <- PC[I8, I8] * outer(kk[I8], kk[I8]); within8 <- within8[upper.tri(within8)]
withinO <- PC[IO, IO][upper.tri(PC[IO, IO])]
# Two-factor model: S + RI on the 8, Obs on its three, correlated (04's syntax with FFMQ Observing alone).
two_syn <- c(sf2b_scale_syntax(I8), sf2b_validity_syntax(OBS, D))
fw2 <- fit_w(two_syn, DAT, c(I8, IO)); f2 <- fw2$fit
stopifnot(!sf2b_failed(f2))
RR <- lavResiduals(f2, type = "cor.bentler")$cov[c(I8, IO), c(I8, IO)]
cross[, residual := RR[cbind(item_8, item_obs)]][, keyed_residual := residual * kk[item_8]]
cross[, rank_keyed_polychoric := frank(-keyed_polychoric)]
setorder(cross, rank_keyed_polychoric)
allres <- as.data.table(as.table(RR)); setnames(allres, c("item_a", "item_b", "residual"))
allres <- allres[match(item_a, c(I8, IO)) < match(item_b, c(I8, IO))]
allres[, block := fifelse(item_a %in% I8 & item_b %in% I8, "within the 8", fifelse(item_a %in% IO & item_b %in% IO, "within Observing", "across"))]
big <- allres[abs(residual) > .10][order(-abs(residual))]
pe2 <- parameterEstimates(f2); r_two <- pe2[pe2$label == "r_ffmq_observing", ]
# Sensitivity: the cross-measure residual covariances above .10 freed, the latent r re-estimated.
big_cross <- big[block == "across"]
sens <- NULL
if (nrow(big_cross)) {
  fw3 <- fit_w(c(two_syn, paste(big_cross$item_a, "~~", big_cross$item_b)), DAT, c(I8, IO))
  if (!sf2b_failed(fw3$fit)) { pe3 <- parameterEstimates(fw3$fit); x <- pe3[pe3$label == "r_ffmq_observing", ]
    sens <- list(freed = paste(big_cross$item_a, "~~", big_cross$item_b), r = x$est, lower = x$ci.lower, upper = x$ci.upper, fit = sf2b_fit_stats(fw3$fit))
  } else sens <- list(freed = paste(big_cross$item_a, "~~", big_cross$item_b), error = sf2b_msg(fw3$fit))
}
secs_c <- as.numeric(difftime(Sys.time(), t_c, units = "secs"))
fwrite(num(copy(cross)), sf2b_out("12_c_item_cross_pairs.csv"))
fwrite(num(copy(big)), sf2b_out("12_c_residuals_above_10.csv"))

## Report ---------
fwrite(fits_a, sf2b_out("12_a_model_fit.csv"))
gv <- function(dr, lab) { x <- share[direction == dr & label == lab]; if (nrow(x)) as.list(x[, .(est, lower, upper)]) else NULL }
out <- list(
  analysis = "12: what is left of the 8 once FFMQ Observing is accounted for. EXPLORATORY, NOT REGISTERED (osf.io/d4a5u does not include it); run after the confirmatory run, with its results seen.",
  mode = SF2B$mode, tag = SF2B$tag, written = format(Sys.time(), "%Y-%m-%d %H:%M"),
  software = list(R = R.version.string, lavaan = as.character(packageVersion("lavaan")), psych = as.character(packageVersion("psych"))),
  n = nrow(DAT), choices = CHOICES,
  a = list(
    r_8_obs_this_model = gv("8|obs", "r_8_obs"), r_8_obs_in_04_ffmq_model = LAT04[score == OBS, est],
    share_of_8_left = gv("8|obs", "unexplained"), share_of_obs_left = gv("obs|8", "unexplained"),
    base_fit = fits_a[model == "base 8|obs"],
    residual_correlations = resid_out, known_groups_latent = mim, known_groups_observed = kg_obs, pause_observed = pause_obs,
    failed = lapply(Filter(function(r) !r$ok, RES), function(r) list(model = paste(r$task$name, r$task$direction), error = r$error))),
  b = lapply(BIF, function(b) b[setdiff(names(b), "loadings")]),
  c = list(cross_pairs = cross, polychoric_within_8_keyed = as.list(summary(within8)), polychoric_within_obs = as.list(summary(withinO)),
           two_factor_fit = sf2b_fit_stats(f2), two_factor_warnings = fw2$warnings,
           r_8_obs_two_factor = list(est = r_two$est, lower = r_two$ci.lower, upper = r_two$ci.upper),
           residuals_above_10 = big, sensitivity_freed_cross_residuals = sens),
  run_seconds = list(a = round(secs_a), b = round(secs_b), c = round(secs_c),
                     total = round(as.numeric(difftime(Sys.time(), T0, units = "secs")))), cores = SF2B$cores)
write_json(out, sf2b_out("12_exploratory_overlap.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
sf2b_log("12 done in ", out$run_seconds$total, " s; outputs ", sf2b_out("12_*"))
print(share); print(resid_out[, .(focal, score, r_total, residual_r, lower, upper)])
print(mim[, .(model, label, est = r3(est), lower = r3(lower), upper = r3(upper))]); print(bif_fit[, .(model, converged, proper, cfi.scaled, rmsea.scaled)])
