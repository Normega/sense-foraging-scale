# 10: decision 46's deviation analysis for the Brief MAIA-2 models. NOT REGISTERED.
#
# This is a deviation from the filed registration (osf.io/d4a5u), fixed by Norm on
# 15 November 2026 (study2/decisions_after_2a.md, decision 46) before any of the values
# it reports had been seen. The registered result stays primary: in it all five Brief
# MAIA-2 validity models fail to converge, so MAIA Noticing is missing, H2 holds at 3 of
# 4 and H3b does not hold (3 of 5). This script is reported beside it, labelled as a
# deviation, and changes no registered cell or output.
#
# The repair, as decision 46 states it. Every model containing maia_12 (the five
# Brief MAIA-2 models of 04_validity.R; the H3d models hold Noticing but not maia_12, so
# there is nothing to change in them) is refitted with one change: maia_12's
# latent-response residual variance is at least .10 of its total variance (standardized
# loading at most .949). Everything else as registered. The first step that works is used:
#   1. the registered settings (WLSMV, theta, pairwise, the agreement factor as
#      registered) with the bound. Under theta, maia_12's residual variance is fixed at 1
#      and Not-Worrying's variance at 1, and maia_12 loads on Not-Worrying only, so its
#      total variance is l^2 + 1 and the bound is l^2 <= 9 (|l| <= 3, std. loading
#      3 / sqrt(10) = .949);
#   2. if lavaan cannot impose it: the delta parameterization in those models only, the
#      same bound (total variance 1, so l^2 <= .90);
#   3. if that also fails: maia_12 left out of Not-Worrying in those models, otherwise
#      as registered.
# A step works when every refitted model converges, is oriented as keyed
# (sf2b_orientation, as 04 requires), and the fitted solution satisfies the bound.
#
# Reported from it, and nothing else: the 8's latent correlation with MAIA Noticing (H2's
# fourth cell, and its upper bound for H3a), contrast 4 (dif_4, as 04 defines it), H3b's
# count with contrast 4 from this repair and contrasts 1, 2, 3 and 5 from the registered
# run, H2's count and H3a with the Noticing cell filled in; fit, convergence, maia_12's
# fitted values and any improper parameter of every refitted model.
#
# Writes only output/<tag>/10_deviation_maia.csv and .json (tag "real", or
# "real_first700" with SF2B_FIRST_N=700, the decision 45 sensitivity analysis). Reads the
# pipeline's saved objects (derived/sf2b_<tag>.rds) and 04's registered decisions.
#
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/10_deviation_maia.R

# Set Up ---------
## Load libraries ---------
packages <- c("jsonlite", "data.table", "lavaan")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

## Settings ---------
source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")
BOUND_SHARE <- .10      # decision 46: residual variance >= .10 of the total
TOL <- 1e-4             # numerical tolerance when checking the fitted bound

# The pipeline's own data and model definitions: 04_validity.R's top-level expressions up
# to (not including) its first fit, which build DAT, TASKS, model_of() and
# contrast_lines() and fit nothing. Its decision helpers (rc_class) are taken the same way.
ex04 <- parse("study2b/analysis/04_validity.R", keep.source = FALSE)
first_fit <- which(vapply(ex04, function(e) grepl("^sf2b_log\\(\"fitting", deparse(e)[1]), TRUE))
stopifnot(length(first_fit) == 1)
for (i in seq_len(first_fit - 1)) eval(ex04[[i]], envir = globalenv())
is_def <- function(e, nm) is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(nm))
for (nm in "rc_class") eval(Filter(function(e) is_def(e, nm), as.list(ex04))[[1]], envir = globalenv())
stopifnot(exists("TASKS"), exists("model_of"), exists("contrast_lines"), exists("DAT"), !exists("RES"))

# The models to refit: every task whose model holds maia_12.
REFIT <- Filter(function(t) "maia_12" %in% model_of(t)$ord, TASKS)
stopifnot(length(REFIT) == 5, all(vapply(REFIT, function(t) t$group == "maia", TRUE)))
task_id <- function(t) sprintf("%s_%s_maia%s", t$form, t$sample, if (t$ri_free) "_rifree" else "")

# The registered estimation call with only the parameterization changed (step 2).
sem_txt <- deparse(sf2b_sem)
stopifnot(sum(grepl("parameterization = \"theta\"", sem_txt, fixed = TRUE)) == 1)
sf2b_sem_delta <- eval(parse(text = sub("parameterization = \"theta\"", "parameterization = \"delta\"", sem_txt, fixed = TRUE)))
environment(sf2b_sem_delta) <- environment(sf2b_sem)

# Step 3's definitions: maia_12 out of Not-Worrying, everything else as built.
D_drop <- D
D_drop$validity$maia$scores$maia_notworrying <- setdiff(D$validity$maia$scores$maia_notworrying, "maia_12")
stopifnot(length(D_drop$validity$maia$scores$maia_notworrying) == 2)

# Model syntax per step: 04's model_of(), plus the bound (steps 1-2) or with the step-3
# definitions (model_of evaluated with D and score_items taken from D_drop).
model_for <- function(t, step) {
  if (step == 3) {
    env <- new.env(parent = environment(model_of)); env$D <- D_drop; env$score_items <- sf2b_score_items(D_drop)
    mf <- model_of; environment(mf) <- env
    return(mf(t))
  }
  m <- model_of(t)
  bound <- if (step == 1) 9 else 1 - BOUND_SHARE    # theta: l^2 <= (1 - .10) / .10 * 1; delta: l^2 <= .90
  m$syn <- c(m$syn, "maia_notworrying =~ l_maia12*maia_12", sprintf("l_maia12^2 < %s", format(bound)))
  m
}

## Fit ---------
improper_of <- function(f) {
  std <- lavInspect(f, "std"); est <- lavInspect(f, "est"); out <- character()
  L <- std$lambda; ix <- which(abs(L) >= 1, arr.ind = TRUE)
  if (nrow(ix)) out <- c(out, sprintf("std loading %s =~ %s = %.3f", colnames(L)[ix[, 2]], rownames(L)[ix[, 1]], L[ix]))
  th <- diag(est$theta); out <- c(out, if (any(th < 0)) sprintf("negative residual variance %s = %.3f", names(th)[th < 0], th[th < 0]))
  sth <- diag(std$theta); out <- c(out, if (any(sth < 0)) sprintf("negative std residual variance %s = %.3f", names(sth)[sth < 0], sth[sth < 0]))
  ps <- diag(est$psi); out <- c(out, if (any(ps < 0)) sprintf("negative factor variance %s = %.3f", names(ps)[ps < 0], ps[ps < 0]))
  C <- lavInspect(f, "cor.lv"); C0 <- C; C0[lower.tri(C0, TRUE)] <- 0
  cx <- which(abs(C0) >= 1, arr.ind = TRUE)
  if (nrow(cx)) out <- c(out, sprintf("factor correlation %s ~~ %s = %.3f", rownames(C)[cx[, 1]], colnames(C)[cx[, 2]], C0[cx]))
  ev <- min(eigen(lavInspect(f, "cov.lv"), symmetric = TRUE, only.values = TRUE)$values)
  if (ev <= 0) out <- c(out, sprintf("factor covariance matrix not positive definite (smallest eigenvalue %.4f)", ev))
  pc <- tryCatch(suppressWarnings(lavInspect(f, "post.check")), error = function(e) NA)
  list(list = out, post_check = isTRUE(pc))
}
fit_one <- function(t, step) {
  m <- model_for(t, step)
  d <- if (t$sample == "unflagged") DAT[!flagged, , drop = FALSE] else DAT
  sem <- if (step == 2) sf2b_sem_delta else sf2b_sem
  f <- sf2b_try(sem(m$syn, d, ordered = m$ord))
  row <- list(model = task_id(t), form = t$form, sample = t$sample, ri_free = t$ri_free, step = step,
              parameterization = if (step == 2) "delta" else "theta", converged = !sf2b_failed(f), error = sf2b_msg(f))
  if (!row$converged) return(c(row, works = FALSE))
  pt <- parTable(f)
  if (step < 3) stopifnot(sum(pt$op == "=~" & pt$rhs == "maia_12") == 1, pt$label[pt$op == "=~" & pt$rhs == "maia_12"] == "l_maia12",
                          sum(pt$op == "<") == 1)
  bad <- sf2b_orientation(f, D[[t$form]], t$scores, if (step == 3) D_drop else D)
  std <- lavInspect(f, "std"); est <- lavInspect(f, "est")
  imp <- improper_of(f)
  fs <- sf2b_fit_stats(f)
  row <- c(row, list(oriented = !length(bad), not_oriented = paste(bad, collapse = " "), optimizer = lavInspect(f, "options")$optim.method,
                     iterations = lavInspect(f, "iterations"),
                     cfi = fs$cfi.scaled, rmsea = fs$rmsea.scaled, rmsea_lower = fs$rmsea.ci.lower.scaled, rmsea_upper = fs$rmsea.ci.upper.scaled,
                     srmr = fs$srmr, chisq = fs$chisq.scaled, df = fs$df, ntotal = fs$ntotal,
                     post_check = imp$post_check, improper = paste(imp$list, collapse = "; ")))
  if (step < 3) {
    l <- est$lambda["maia_12", "maia_notworrying"]; sl <- std$lambda["maia_12", "maia_notworrying"]
    share <- diag(std$theta)[["maia_12"]]
    row <- c(row, list(maia12_loading = l, maia12_std_loading = sl, maia12_residual_variance = diag(est$theta)[["maia_12"]],
                       maia12_std_residual_variance = share, bound_held = share >= BOUND_SHARE - TOL, bound_active = share < BOUND_SHARE + TOL))
  }
  row$works <- row$converged && row$oriented && (step == 3 || isTRUE(row$bound_held))
  pe <- parameterEstimates(f)
  keep <- pe$label %in% c("r_maia_noticing", "r_maia_notworrying", "dif_4")
  row$estimates <- pe[keep, c("label", "est", "se", "ci.lower", "ci.upper")]
  row
}

RUN <- list(); STEP <- NA_integer_
for (step in 1:3) {
  sf2b_log("decision 46, step ", step, ": refitting ", length(REFIT), " Brief MAIA-2 models")
  res <- sf2b_par_lapply(REFIT, function(t) fit_one(t, STEP_NOW),
                         data = list(STEP_NOW = step, fit_one = fit_one, model_for = model_for, model_of = model_of, improper_of = improper_of,
                                     task_id = task_id, sf2b_sem_delta = sf2b_sem_delta, DAT = DAT, D = D, D_drop = D_drop,
                                     score_items = score_items, flagged = M$flagged, contrast_lines = contrast_lines,
                                     BOUND_SHARE = BOUND_SHARE, TOL = TOL))
  RUN[[step]] <- res
  for (r in res) sf2b_log("  ", r$model, ": converged ", r$converged, if (nzchar(r$error)) paste0(" (", r$error, ")"),
                          if (!is.null(r$maia12_std_loading)) sprintf("; maia_12 std loading %.4f, residual share %.4f", r$maia12_std_loading, r$maia12_std_residual_variance),
                          if (!is.null(r$works)) paste0("; works ", r$works))
  if (all(vapply(res, function(r) isTRUE(r$works), TRUE))) { STEP <- step; break }
}

## Report ---------
fit_tab <- rbindlist(lapply(unlist(RUN, recursive = FALSE), function(r) as.data.table(r[setdiff(names(r), "estimates")])), fill = TRUE)
fit_tab[, used := !is.na(STEP) & step == STEP]
num <- names(fit_tab)[vapply(fit_tab, is.numeric, TRUE) & !names(fit_tab) %in% c("step", "iterations", "df", "ntotal")]
fit_tab[, (num) := lapply(.SD, function(x) round(x, 4)), .SDcols = num]

reg <- fromJSON(sf2b_out("04_validity_decisions.json"), simplifyVector = TRUE)$sf8
out <- list(analysis = "decision 46 deviation analysis (NOT registered): Brief MAIA-2 models refitted with maia_12's residual variance >= .10 of its total",
            mode = SF2B$mode, tag = SF2B$tag, written = format(Sys.time(), "%Y-%m-%d %H:%M"),
            software = list(R = R.version.string, lavaan = as.character(packageVersion("lavaan"))),
            step_used = STEP,
            step_note = c("1" = "registered settings (theta) with the bound", "2" = "delta parameterization with the bound",
                          "3" = "maia_12 left out of Not-Worrying")[as.character(STEP)],
            steps_tried = lapply(RUN, function(res) lapply(res, function(r) r[c("model", "converged", "error", "works")])))
if (!is.na(STEP)) {
  prim <- Filter(function(r) r$model == "sf8_all_maia", RUN[[STEP]])[[1]]
  e <- as.data.table(prim$estimates)
  g <- function(lab) e[label == lab]
  nt <- g("r_maia_noticing"); d4 <- g("dif_4"); nw <- g("r_maia_notworrying")
  h2 <- as.data.table(reg$h2)
  h2_cell <- list(score = "maia_noticing", r = r3(nt$est), lower = r3(nt$ci.lower), upper = r3(nt$ci.upper), se = r3(nt$se),
                  h2 = isTRUE(nt$ci.lower > SF2B$h2_lower), h3a = isTRUE(nt$ci.upper < SF2B$h3a_upper),
                  rc_class = rc_class(nt$est, nt$ci.lower, nt$ci.upper))
  h2_other <- h2[score != "maia_noticing"]
  h2_n <- sum(h2_other$h2) + h2_cell$h2
  h3b <- as.data.table(reg$h3b)
  c4 <- list(contrast = 4L, first = "maia_noticing", second = "maia_notworrying", r_first = r3(nt$est), r_second = r3(nw$est),
             difference = r3(d4$est), se = r3(d4$se), lower = r3(d4$ci.lower), upper = r3(d4$ci.upper), holds = isTRUE(d4$ci.lower > 0))
  h3b_n <- sum(h3b[contrast != 4]$holds) + c4$holds
  out$sf8 <- list(
    h2_noticing = h2_cell,
    H2 = list(correlations_holding = h2_n, of = 4L, needed = SF2B$h2_needed, holds = h2_n >= SF2B$h2_needed,
              registered_cells = h2_other[, .(score, r, lower, upper, h2)]),
    H3a = list(holds = all(h2_other$h3a) && h2_cell$h3a, upper_bounds = as.list(c(setNames(h2_other$upper, h2_other$score), maia_noticing = h2_cell$upper)),
               rc_classes = as.list(c(setNames(h2_other$rc_class, h2_other$score), maia_noticing = h2_cell$rc_class))),
    contrast_4 = c4,
    H3b = list(contrasts_holding = h3b_n, of = 5L, needed = SF2B$h3b_needed, holds = h3b_n >= SF2B$h3b_needed,
               registered_contrasts = h3b[contrast != 4, .(contrast, first, second, difference, lower, upper, holds)]),
    maia_12 = prim[c("maia12_loading", "maia12_std_loading", "maia12_residual_variance", "maia12_std_residual_variance", "bound_held", "bound_active")])
}
out$models <- fit_tab[used == TRUE]
fwrite(fit_tab, sf2b_out("10_deviation_maia.csv"))
write_json(out, sf2b_out("10_deviation_maia.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
sf2b_log("step used: ", STEP, "; written ", sf2b_out("10_deviation_maia.csv"), " and .json")
if (!is.na(STEP)) print(out$sf8[c("h2_noticing", "contrast_4")])
