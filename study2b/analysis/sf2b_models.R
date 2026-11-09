# Model helpers shared by the Study 2b scripts, so every refit (subsets, bootstrap
# resamples, the 10) is the identical model.
#
# Registered settings in every fit (Section 5, "Statistical models"): items ordinal,
# WLSMV, theta parameterization, missing = "pairwise" (set explicitly); scaled fit
# indices; the agreement factor (RI) a random intercept with loadings fixed at 1 on
# every UNRECODED statement of the scale in the model, free variance, uncorrelated with
# the content factor and with every validity factor (decision 5).
# conditional.x = FALSE everywhere (Study 2's finding: with ordinal indicators lavaan's
# default gave wrong regression estimates whenever an observed covariate is present).

# ---- syntax ---------------------------------------------------------------------------------
# The H1 model: one content factor S (variance 1, loadings free) and the RI.
sf2b_scale_syntax <- function(items, ri = TRUE) {
  s <- c(paste0("S =~ NA*", paste(items, collapse = " + ")), "S ~~ 1*S")
  if (ri) s <- c(s, paste0("RI =~ ", paste0("1*", items, collapse = " + ")), "RI ~~ RI", "S ~~ 0*RI")
  s
}
# Validity factors: one per score in `scores`, each on its items (unrecoded: a reverse-
# keyed validity item simply loads negatively), variance 1, correlated with S (labelled
# r_<score>) and with each other, uncorrelated with the RI unless ri_free.
# A factor with free loadings and variance 1 has an arbitrary sign; lavaan's starts
# orient it by its first item, which for BFI-2 Aesthetic Sensitivity, SBI and several
# FFMQ and Brief MAIA-2 subscales is reverse-keyed (found in the simulated dry run:
# their correlations came out with the wrong sign). So every loading starts at +/-.6 by
# its keying, and sf2b_orientation() checks the result.
sf2b_validity_syntax <- function(scores, defs, ri_free = FALSE, ri = TRUE) {
  its <- sf2b_score_items(defs); rv <- unlist(lapply(defs$validity, `[[`, "reverse"))
  unlist(lapply(scores, function(sc) c(
    paste0(sc, " =~ NA*", its[[sc]][1], " + ", paste0("start(", ifelse(its[[sc]] %in% rv, -.6, .6), ")*", its[[sc]], collapse = " + ")),
    paste0(sc, " ~~ 1*", sc),
    paste0("S ~~ r_", sc, "*", sc),
    if (ri) paste0("RI ~~ ", if (ri_free) paste0("ri_", sc, "*") else "0*", sc))))
}
# Every factor oriented as keyed: the sum of its keyed standardized loadings is positive
# (S over the scale's statements, each validity factor over its items). Returns the
# names of factors that are not.
sf2b_orientation <- function(fit, items, scores, defs) {
  L <- lavInspect(fit, "std")$lambda; its <- sf2b_score_items(defs)
  rv <- c(defs$reversed, unlist(lapply(defs$validity, `[[`, "reverse")))
  keyed_sum <- function(f, it) sum(L[it, f] * ifelse(it %in% rv, -1, 1))
  bad <- if (keyed_sum("S", items) <= 0) "S" else character()
  for (sc in scores) if (keyed_sum(sc, its[[sc]]) <= 0) bad <- c(bad, sc)
  bad
}
# score name -> its items, over every validity measure.
sf2b_score_items <- function(defs) do.call(c, unname(lapply(defs$validity, `[[`, "scores")))

sf2b_sem <- function(model, data, ordered, fast = FALSE, group = NULL, ...) {
  args <- list(model = paste(model, collapse = "\n"), data = data, ordered = ordered, estimator = "WLSMV",
               parameterization = "theta", missing = "pairwise", conditional.x = FALSE, ...)
  if (fast) args <- c(args, list(se = "none", test = "none"))
  if (!is.null(group)) args$group <- group
  suppressWarnings(do.call(lavaan::sem, args))
}
sf2b_try <- function(expr) tryCatch(expr, error = function(e) e)
sf2b_failed <- function(x) inherits(x, "error") || (inherits(x, "lavaan") && !lavInspect(x, "converged"))
sf2b_msg <- function(x) if (inherits(x, "error")) trimws(gsub("[[:space:]]+", " ", conditionMessage(x))) else if (sf2b_failed(x)) "did not converge" else ""

# ---- H1 (decisions 17, 33) ----------------------------------------------------------------------
# Content omega from the standardized solution, keyed: (sum l)^2 / ((sum l)^2 + (sum r)^2
# + sum(Theta)), l = keyed content loadings, r = keyed agreement loadings, Theta the
# keyed residual (co)variance matrix (so a freed residual covariance counts as error).
# The agreement factor counts as error. This is the formula of the power simulation
# (study2/power/power_2b.R), which the registered rule was calibrated on.
sf2b_omega <- function(fit, items, reversed) {
  s <- lavInspect(fit, "std"); k <- ifelse(items %in% reversed, -1, 1)
  l <- s$lambda[items, "S"] * k
  r <- if ("RI" %in% colnames(s$lambda)) s$lambda[items, "RI"] * k else 0
  Th <- s$theta[items, items] * outer(k, k)
  sum(l)^2 / (sum(l)^2 + sum(r)^2 + sum(Th))
}

SF2B_FIT_INDICES <- c("ntotal", "npar", "chisq.scaled", "df", "pvalue.scaled", "cfi.scaled", "tli.scaled", "rmsea.scaled",
                      "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "srmr", "cfi.robust", "rmsea.robust")
sf2b_fit_stats <- function(fit) {
  if (sf2b_failed(fit)) return(setNames(as.list(rep(NA_real_, length(SF2B_FIT_INDICES))), SF2B_FIT_INDICES))
  as.list(round(unclass(suppressWarnings(fitMeasures(fit, SF2B_FIT_INDICES))), 4))
}

# The registered rule on one item set. Returns every quantity the rule reads, and the decision.
# PREREG: "the largest" residual is the largest in absolute value; freeing it adds that
# residual covariance to the same model; "every standardized content loading" is
# compared item by item between the two fits (std.all, keyed sign kept as is).
sf2b_h1 <- function(X, items, reversed) {
  syn <- sf2b_scale_syntax(items)
  f <- sf2b_try(sf2b_sem(syn, as.data.frame(X[, items]), ordered = items))
  out <- list(converged = !sf2b_failed(f), error = sf2b_msg(f))
  if (!out$converged) return(c(out, list(supported = FALSE)))
  fs <- sf2b_fit_stats(f)
  rr <- lavResiduals(f, type = "cor.bentler", summary = TRUE)
  R <- rr$cov[items, items]; R0 <- R; R0[upper.tri(R0, TRUE)] <- NA
  ij <- arrayInd(which.max(abs(R0)), dim(R0))
  pair <- c(items[ij[1]], items[ij[2]])
  maxres <- max(abs(R0), na.rm = TRUE)
  om <- sf2b_omega(f, items, reversed)
  f2 <- sf2b_try(sf2b_sem(c(syn, paste(pair[1], "~~", pair[2])), as.data.frame(X[, items]), ordered = items))
  d_om <- d_lam <- NA_real_
  if (!sf2b_failed(f2)) {
    d_om <- abs(sf2b_omega(f2, items, reversed) - om)
    d_lam <- max(abs(lavInspect(f2, "std")$lambda[items, "S"] - lavInspect(f, "std")$lambda[items, "S"]))
  }
  sm <- rr$summary
  closefit_p <- sm["usrmr.closefit.pvalue", "cor"]
  cfi_ok <- isTRUE(fs$cfi.scaled >= SF2B$h1_cfi)
  local_ok <- maxres <= SF2B$h1_resid || isTRUE(d_om < SF2B$h1_d_omega && d_lam < SF2B$h1_d_lambda)
  close_ok <- isTRUE(closefit_p > SF2B$h1_closefit_alpha)
  std <- lavInspect(f, "std")
  pe <- parameterEstimates(f)
  riv <- pe[pe$lhs == "RI" & pe$op == "~~" & pe$rhs == "RI", ]
  c(out, list(fit = f, stats = fs, max_residual = maxres, max_residual_pair = paste(pair, collapse = " ~~ "),
              d_omega = d_om, d_lambda = d_lam, usrmr = sm["usrmr", "cor"], usrmr_ci = unname(sm[c("usrmr.ci.lower", "usrmr.ci.upper"), "cor"]),
              closefit_p = closefit_p, omega = om, cfi_ok = cfi_ok, local_ok = local_ok, close_ok = close_ok,
              supported = cfi_ok && local_ok && close_ok, residuals = R,
              loadings = data.table(item = items, content = std$lambda[items, "S"], agreement = std$lambda[items, "RI"],
                                    keyed_content = std$lambda[items, "S"] * ifelse(items %in% reversed, -1, 1),
                                    residual_variance = diag(std$theta)[items]),
              ri_variance = riv$est, ri_se = riv$se, ri_p = riv$pvalue,
              ri_share = mean(std$lambda[items, "RI"]^2)))
}

# ---- observed scores --------------------------------------------------------------------------------
# Mean of answered items after reverse-scoring ((min + max) - x) when at least `need`
# are answered (default: the 80% rule); "sum" scales are the mean x k (a prorated sum).
sf2b_score <- function(X, items, reversed = character(), range = c(1, 6), method = "mean", need = NULL) {
  M <- X[, items, drop = FALSE]; storage.mode(M) <- "double"
  rv <- intersect(items, reversed)
  if (length(rv)) M[, rv] <- (range[1] + range[2]) - M[, rv]
  if (is.null(need)) need <- ceiling(SF2B$min_prop_answered * length(items) - 1e-9)
  s <- rowMeans(M, na.rm = TRUE); s[rowSums(!is.na(M)) < need] <- NA
  if (method == "sum") s <- s * length(items)
  s
}

# ---- known groups -------------------------------------------------------------------------------------
# lm(score ~ group + covariates); d = the group coefficient over the pooled SD of the
# score within the two groups; its 95% CI is the coefficient's CI over that SD. With no
# covariates the coefficient's t-test is Student's two-sample t-test.
sf2b_kg <- function(y, g, covars = NULL) {
  d <- data.table(y = y, g = as.integer(g)); if (!is.null(covars)) d <- cbind(d, covars)
  d <- d[complete.cases(d)]
  cv <- setdiff(names(d), c("y", "g"))
  cv <- cv[vapply(cv, function(v) length(unique(d[[v]])) > 1, TRUE)]
  m <- lm(reformulate(c("g", cv), "y"), data = d)
  s <- summary(m)$coefficients["g", ]; ci <- confint(m)["g", ]
  n1 <- sum(d$g == 1); n0 <- sum(d$g == 0)
  sp <- sqrt(((n1 - 1) * var(d$y[d$g == 1]) + (n0 - 1) * var(d$y[d$g == 0])) / (n1 + n0 - 2))
  data.table(n_group = n1, n_reference = n0, mean_group = mean(d$y[d$g == 1]), mean_reference = mean(d$y[d$g == 0]),
             coef = s[["Estimate"]], se = s[["Std. Error"]], t = s[["t value"]], p = s[["Pr(>|t|)"]], pooled_sd = sp,
             d = s[["Estimate"]] / sp, d_lower = ci[[1]] / sp, d_upper = ci[[2]] / sp,
             covariates = paste(cv, collapse = " + "))
}

sf2b_r_ci <- function(x, y, method = "pearson") {
  ok <- complete.cases(x, y)
  if (sum(ok) < 4 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(list(n = sum(ok), r = NA_real_, lower = NA_real_, upper = NA_real_, p = NA_real_))
  t <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
  ci <- if (!is.null(t$conf.int)) t$conf.int else {   # Spearman: Fisher z with the Bonett-Wright SE
    z <- atanh(t$estimate); se <- sqrt((1 + t$estimate^2 / 2) / (sum(ok) - 3)); tanh(z + c(-1, 1) * qnorm(.975) * se) }
  list(n = sum(ok), r = unname(t$estimate), lower = ci[1], upper = ci[2], p = t$p.value)
}

# ---- parallel map and bootstrap ----------------------------------------------------------------------
# FUN runs on PSOCK workers that have config.R and this file sourced; `data` becomes
# FUN's environment (and that of any function in `data`), so nothing else from the
# calling script is serialized.
sf2b_par_lapply <- function(X, FUN, data = list(), cores = SF2B$cores) {
  env <- list2env(data, parent = globalenv())
  for (nm in names(data)) if (is.function(data[[nm]])) { f <- data[[nm]]; environment(f) <- env; assign(nm, f, envir = env) }
  environment(FUN) <- env
  if (cores <= 1 || length(X) < 2) return(lapply(X, FUN))
  cl <- makeCluster(min(cores, length(X)))
  on.exit(stopCluster(cl))
  clusterCall(cl, function(wd, env) { setwd(wd); do.call(Sys.setenv, as.list(env)); NULL }, getwd(),
              c(SF2B_MODE = SF2B$mode, SF2B_DATA_DIR = SF2B$data_dir, SF2B_CORES = "1"))
  clusterEvalQ(cl, { suppressMessages({ source("study2b/analysis/config.R"); source("study2b/analysis/sf2b_models.R") }); NULL })
  parLapplyLB(cl, X, FUN)
}
# Index sets drawn in the calling process with a recorded seed, so results do not
# depend on how workers are scheduled.
sf2b_boot_indices <- function(n, B, seed) { set.seed(seed); lapply(seq_len(B), function(b) sample.int(n, n, replace = TRUE)) }

# ---- invariance (exploratory: H8, gender, country; Wu & Estabrook 2016, theta) -----------------------
# Sparse categories: a response category with fewer than 5 answers in ANY group is merged
# with its neighbour toward the scale's centre, for that item, in all groups, until none
# is sparse (or two remain); categories are then renumbered 1..m.
sf2b_merge_sparse <- function(X, group, min_n = SF2B$min_category_n, centre = 3.5) {
  group <- factor(group); log <- list()
  for (j in colnames(X)) {
    x <- X[, j]
    repeat {
      cats <- sort(unique(x[!is.na(x)]))
      if (length(cats) <= 2) break
      tab <- table(factor(x, levels = cats), group)
      sparse <- cats[apply(tab, 1, min) < min_n]
      if (!length(sparse)) break
      c0 <- sparse[which.max(abs(sparse - centre))]
      into <- if (c0 < centre) min(cats[cats > c0]) else max(cats[cats < c0])
      log[[length(log) + 1]] <- data.table(item = j, category = c0, merged_into = into)
      x[!is.na(x) & x == c0] <- into
    }
    X[, j] <- match(x, sort(unique(x[!is.na(x)])))
  }
  list(X = X, log = rbindlist(log))
}

# One content factor S plus the RI, G groups. Steps (Section 5, H8): configural
# (thresholds group-specific; latent-response intercepts 0, residual variances 1; factor
# mean 0, variance 1), thresholds (equal; intercepts and residual variances free in
# groups 2..G), loadings (+ equal; S's variance free in groups 2..G). The RI keeps
# loadings 1 in every group, its variance free in every group, and is uncorrelated with
# S. `free`: items released from a step's constraint (partial invariance).
SF2B_INV_STEPS <- c("configural", "thresholds", "loadings")
sf2b_inv_syntax <- function(items, ncat, G, step, free = list()) {
  lev <- match(step, SF2B_INV_STEPS)
  gl <- function(x) paste0("c(", paste(x, collapse = ", "), ")")
  grp <- function(lab) paste0(lab, "_g", seq_len(G))
  terms <- vapply(items, function(i) {
    eq <- lev >= 3 && !(i %in% free$loadings)
    paste0(gl(if (eq) rep(paste0("l_", i), G) else grp(paste0("l_", i))), "*", i) }, "")
  s <- c(paste0("S =~ ", paste(terms, collapse = " + ")), paste0("RI =~ ", paste0("1*", items, collapse = " + ")))
  for (i in items) {
    thr_eq <- lev >= 2 && !(i %in% free$thresholds)
    for (k in seq_len(ncat[[i]] - 1))
      s <- c(s, paste0(i, " | ", gl(if (thr_eq) rep(paste0("t_", i, "_", k), G) else grp(paste0("t_", i, "_", k))), "*t", k))
    can_free <- thr_eq && ncat[[i]] >= 3
    s <- c(s, paste0(i, " ~ ", gl(if (can_free) c("0", rep("NA", G - 1)) else rep("0", G)), "*1"),
           paste0(i, " ~~ ", gl(if (can_free) c("1", rep("NA", G - 1)) else rep("1", G)), "*", i))
  }
  s <- c(s, paste0("S ~~ ", gl(if (lev >= 3) c("1", rep("NA", G - 1)) else rep("1", G)), "*S"),
         paste0("S ~ ", gl(rep("0", G)), "*1"),
         paste0("RI ~~ ", gl(rep("NA", G)), "*RI"), paste0("RI ~ ", gl(rep("0", G)), "*1"),
         paste0("RI ~~ ", gl(rep("0", G)), "*S"))
  paste(s, collapse = "\n")
}

# Runs the steps in order; a step is non-invariant when delta scaled CFI < -.010 AND
# delta scaled RMSEA > .015 (PREREG: "and" as written in Section 5, H8). Then partial
# invariance: each item still constrained at that step is released in turn, the release
# with the best CFI kept, until the step passes or SF2B$partial_max_free are released.
sf2b_run_mg <- function(X, items, group, label) {
  keep <- !is.na(group)
  mg <- sf2b_merge_sparse(X[keep, items, drop = FALSE], group[keep])
  ncat <- apply(mg$X, 2, function(x) length(unique(x[!is.na(x)])))
  d <- as.data.frame(mg$X); d$.group <- as.character(group[keep]); G <- length(unique(d$.group))
  d <- d[order(d$.group), ]   # lavaan takes groups in order of appearance: the alphabetically first is the reference
  fit_one <- function(step, free, start = NULL) {
    syn <- sf2b_inv_syntax(items, ncat, G, step, free)
    a <- list(syn, d, items, group = ".group", auto.fix.first = FALSE, control = list(iter.max = 3000))
    if (inherits(start, "lavaan")) a$start <- start
    f <- sf2b_try(do.call(sf2b_sem, a))
    if (sf2b_failed(f) && !is.null(a$start)) { a$start <- NULL; f <- sf2b_try(do.call(sf2b_sem, a)) }
    f
  }
  row <- function(f, st, released, prev) {
    fs <- sf2b_fit_stats(f)
    data.table(analysis = label, step = st, released = paste(released, collapse = " "), converged = !sf2b_failed(f),
               error = sf2b_msg(f), chisq = fs$chisq.scaled, df = fs$df, cfi = fs$cfi.scaled, rmsea = fs$rmsea.scaled,
               srmr = fs$srmr, d_cfi = if (!is.null(prev)) round(fs$cfi.scaled - prev$cfi, 4) else NA_real_,
               d_rmsea = if (!is.null(prev)) round(fs$rmsea.scaled - prev$rmsea, 4) else NA_real_)
  }
  noninv <- function(r) !r$converged | (!is.na(r$d_cfi) & r$d_cfi < SF2B$inv_d_cfi & r$d_rmsea > SF2B$inv_d_rmsea)
  sf2b_log("invariance: ", label, " (", paste(names(table(d$.group)), table(d$.group), sep = " = ", collapse = ", "), ")")
  res <- list(); free <- list(thresholds = character(), loadings = character()); prev <- NULL; prev_fit <- NULL
  for (st in SF2B_INV_STEPS) {
    f <- fit_one(st, free, prev_fit); r <- row(f, st, character(), prev)
    if (st == "configural") { res[[1]] <- r; if (!r$converged) break; prev_fit <- f
      fs <- sf2b_fit_stats(f); prev <- list(cfi = fs$cfi.scaled, rmsea = fs$rmsea.scaled); next }
    r[, noninvariant := noninv(r)]; res[[length(res) + 1]] <- r
    if (isTRUE(r$noninvariant)) {
      released <- character()
      repeat {
        if (length(released) >= SF2B$partial_max_free) break
        cand <- setdiff(items, c(released, free[[st]]))
        cf <- vapply(cand, function(i) { fr <- free; fr[[st]] <- c(fr[[st]], released, i)
          f2 <- fit_one(st, fr, prev_fit); if (sf2b_failed(f2)) NA_real_ else fitMeasures(f2, "cfi.scaled")[[1]] }, 0)
        if (all(is.na(cf))) break
        released <- c(released, cand[which.max(cf)])
        fr <- free; fr[[st]] <- c(fr[[st]], released); f <- fit_one(st, fr, prev_fit)
        r <- row(f, paste0(st, " (partial)"), released, prev); r[, noninvariant := noninv(r)]
        res[[length(res) + 1]] <- r
        if (!isTRUE(r$noninvariant)) break
      }
      free[[st]] <- c(free[[st]], released)
      if (isTRUE(r$noninvariant)) break
    }
    prev_fit <- f; fs <- sf2b_fit_stats(f); prev <- list(cfi = fs$cfi.scaled, rmsea = fs$rmsea.scaled)
  }
  list(table = rbindlist(res, fill = TRUE), free = free, groups = as.list(table(d$.group)),
       merges = if (nrow(mg$log)) cbind(analysis = label, mg$log) else NULL)
}
